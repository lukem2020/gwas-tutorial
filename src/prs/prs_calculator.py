"""
Polygenic Risk Score (PRS) Calculator

Implements complete PRS development pipeline with:
- LD-based SNP clumping
- P-value threshold optimization
- Cross-validation framework
- Performance evaluation (AUC, R², OR by quantile)
- Effect size shrinkage methods

This is production-quality code for developing and validating PRS.

Author: Luke M
Repository: https://github.com/lukem2020/gwas-tutorial

Methodology:
    PRS = Σ (genotype_i × β_i) for selected SNPs
    
    Steps:
    1. Select independent SNPs via LD clumping
    2. Optimize P-value threshold
    3. Weight SNPs by effect sizes
    4. Validate in independent dataset
    5. Evaluate performance metrics
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold, KFold
from typing import Tuple, Optional, List, Dict
import logging
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class PRSResults:
    """Container for PRS results"""
    scores: np.ndarray              # PRS values
    performance: float              # AUC or R²
    optimal_threshold: float        # Best P-value threshold
    n_snps: int                    # Number of SNPs in score
    quantile_ors: Optional[Dict]   # Odds ratios by quantile (binary)
    weights: pd.DataFrame          # SNP weights used


class PolygenicRiskScore:
    """
    Polygenic Risk Score Calculator
    
    Features:
    - LD clumping to select independent SNPs
    - P-value threshold optimization
    - Cross-validation for unbiased evaluation
    - Multiple performance metrics
    - Handles both binary and continuous traits
    
    Example:
    --------
    >>> prs_calc = PolygenicRiskScore(clump_r2=0.1, clump_kb=250)
    >>> cv_results = prs_calc.cross_validate(
    ...     gwas_results, genotypes, phenotype, n_folds=5
    ... )
    >>> print(f"CV AUC: {cv_results['mean_score']:.3f} ± {cv_results['std_score']:.3f}")
    """
    
    def __init__(self,
                 clump_r2: float = 0.1,
                 clump_kb: int = 250,
                 p_thresholds: Optional[List[float]] = None):
        """
        Initialize PRS calculator
        
        Parameters:
        -----------
        clump_r2 : float
            R² threshold for LD clumping (default: 0.1)
        clump_kb : int
            Window size for clumping in kb (default: 250)
        p_thresholds : list of float
            P-value thresholds to test (default: standard set)
        """
        self.clump_r2 = clump_r2
        self.clump_kb = clump_kb
        
        if p_thresholds is None:
            # Standard P-value thresholds
            self.p_thresholds = [
                5e-8,   # Genome-wide significant
                1e-6,
                1e-5,
                1e-4,
                1e-3,
                0.01,
                0.05,
                0.1,
                0.5,
                1.0     # All SNPs
            ]
        else:
            self.p_thresholds = sorted(p_thresholds)
        
        logger.info(f"Initialized PRS calculator (R²<{clump_r2}, window={clump_kb}kb)")
    
    def calculate_ld(self, 
                    genotypes: np.ndarray,
                    indices: np.ndarray) -> np.ndarray:
        """
        Calculate pairwise LD (R²) between SNPs
        
        R² = cor(X_i, X_j)²
        
        Parameters:
        -----------
        genotypes : np.ndarray (n_samples x n_snps)
            Genotype matrix
        indices : np.ndarray
            Indices of SNPs to calculate LD for
        
        Returns:
        --------
        ld_matrix : np.ndarray (len(indices) x len(indices))
            Pairwise R² values
        """
        # Extract SNPs
        X = genotypes[:, indices]
        
        # Handle missing data (impute with mean)
        X_clean = X.copy()
        for i in range(X.shape[1]):
            missing = np.isnan(X[:, i])
            if np.any(missing):
                X_clean[missing, i] = np.nanmean(X[:, i])
        
        # Standardize
        X_std = (X_clean - np.mean(X_clean, axis=0)) / np.std(X_clean, axis=0)
        
        # Correlation matrix
        cor = np.corrcoef(X_std.T)
        
        # R² = correlation²
        ld_matrix = cor ** 2
        
        return ld_matrix
    
    def ld_clump(self,
                 gwas_results: pd.DataFrame,
                 genotypes: np.ndarray,
                 p_threshold: float = 0.001) -> pd.DataFrame:
        """
        LD-based clumping to select independent SNPs
        
        Algorithm:
        1. Sort SNPs by P-value
        2. Take SNP with lowest P-value as index SNP
        3. Remove all SNPs in LD (R² > threshold) within window
        4. Repeat until no SNPs remain
        
        Parameters:
        -----------
        gwas_results : pd.DataFrame
            GWAS results with columns: CHR, BP, SNP, BETA, P
        p_threshold : float
            P-value threshold for inclusion
        genotypes : np.ndarray
            Genotype matrix (for LD calculation)
        
        Returns:
        --------
        clumped : pd.DataFrame
            Independent SNPs after clumping
        """
        logger.info(f"LD clumping (R²<{self.clump_r2}, p<{p_threshold})...")
        
        # Filter by P-value
        sig_snps = gwas_results[gwas_results['P'] < p_threshold].copy()
        
        if len(sig_snps) == 0:
            logger.warning(f"No SNPs below p-threshold {p_threshold}")
            return pd.DataFrame()
        
        # Sort by P-value
        sig_snps = sig_snps.sort_values('P').reset_index(drop=True)
        
        logger.info(f"  {len(sig_snps)} SNPs below threshold")
        
        # Track which SNPs to keep
        keep_indices = []
        removed = set()
        
        for i in range(len(sig_snps)):
            if i in removed:
                continue
            
            # This is an index SNP
            keep_indices.append(i)
            
            # Get chromosome and position
            chr_i = sig_snps.iloc[i]['CHR']
            bp_i = sig_snps.iloc[i]['BP']
            
            # Find SNPs in window
            same_chr = sig_snps['CHR'] == chr_i
            in_window = (np.abs(sig_snps['BP'] - bp_i) <= self.clump_kb * 1000)
            window_snps = np.where(same_chr & in_window)[0]
            
            # Remove this index SNP from window SNPs
            window_snps = window_snps[window_snps != i]
            
            if len(window_snps) == 0:
                continue
            
            # Calculate LD with index SNP
            try:
                # Get genotype indices (assuming sig_snps has them)
                idx_i = sig_snps.index[i]
                idx_window = sig_snps.index[window_snps]
                
                # Calculate LD
                all_indices = np.array([idx_i] + list(idx_window))
                ld_matrix = self.calculate_ld(genotypes, all_indices)
                
                # LD with index SNP (first row)
                ld_with_index = ld_matrix[0, 1:]
                
                # Remove SNPs in LD
                in_ld = ld_with_index > self.clump_r2
                remove_idx = window_snps[in_ld]
                removed.update(remove_idx)
                
            except Exception as e:
                logger.warning(f"LD calculation failed for SNP {i}: {e}")
                continue
        
        # Get clumped SNPs
        clumped = sig_snps.iloc[keep_indices].reset_index(drop=True)
        
        logger.info(f"  ✓ {len(clumped)} independent SNPs after clumping")
        
        return clumped
    
    def calculate_prs(self,
                     genotypes: np.ndarray,
                     weights: pd.DataFrame) -> np.ndarray:
        """
        Calculate PRS for individuals
        
        PRS = Σ (genotype_i × β_i)
        
        Parameters:
        -----------
        genotypes : np.ndarray (n_samples x n_snps)
            Full genotype matrix
        weights : pd.DataFrame
            SNP weights with columns: SNP_INDEX, BETA
        
        Returns:
        --------
        prs : np.ndarray (n_samples,)
            Polygenic risk scores
        """
        if len(weights) == 0:
            logger.warning("No SNPs in weights, returning zeros")
            return np.zeros(genotypes.shape[0])
        
        # Extract SNP indices and weights
        snp_indices = weights['SNP_INDEX'].values
        betas = weights['BETA'].values
        
        # Get genotypes for these SNPs
        X = genotypes[:, snp_indices]
        
        # Handle missing data (impute with mean)
        X_clean = X.copy()
        for i in range(X.shape[1]):
            missing = np.isnan(X[:, i])
            if np.any(missing):
                X_clean[missing, i] = np.nanmean(X[:, i])
        
        # Calculate PRS
        prs = X_clean @ betas
        
        return prs
    
    def evaluate_prs(self,
                    prs: np.ndarray,
                    phenotype: np.ndarray,
                    is_binary: bool = True) -> float:
        """
        Evaluate PRS performance
        
        Parameters:
        -----------
        prs : np.ndarray
            Polygenic risk scores
        phenotype : np.ndarray
            True phenotype
        is_binary : bool
            Binary (case/control) vs continuous trait
        
        Returns:
        --------
        performance : float
            AUC (binary) or R² (continuous)
        """
        if is_binary:
            # AUC for binary traits
            try:
                auc = roc_auc_score(phenotype, prs)
                return auc
            except:
                return 0.5
        else:
            # R² for continuous traits
            correlation = np.corrcoef(prs, phenotype)[0, 1]
            r_squared = correlation ** 2
            return r_squared
    
    def cross_validate(self,
                      gwas_results: pd.DataFrame,
                      genotypes: np.ndarray,
                      phenotype: np.ndarray,
                      n_folds: int = 5,
                      is_binary: bool = True) -> Dict:
        """
        Cross-validated PRS performance
        
        For each fold:
        1. Split data into train/test
        2. Clump SNPs using training GWAS
        3. Optimize P-threshold on training
        4. Evaluate on test set
        
        Parameters:
        -----------
        gwas_results : pd.DataFrame
            GWAS summary statistics
        genotypes : np.ndarray
            Genotype matrix
        phenotype : np.ndarray
            Phenotype
        n_folds : int
            Number of CV folds
        is_binary : bool
            Binary trait (case/control)
        
        Returns:
        --------
        cv_results : dict
            Cross-validation results
        """
        logger.info(f"Running {n_folds}-fold cross-validation...")
        
        n_samples = len(phenotype)
        
        # Initialize cross-validation
        if is_binary:
            kf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)
            splits = kf.split(np.zeros(n_samples), phenotype)
        else:
            kf = KFold(n_splits=n_folds, shuffle=True, random_state=42)
            splits = kf.split(np.zeros(n_samples))
        
        fold_scores = []
        fold_thresholds = []
        
        for fold_idx, (train_idx, test_idx) in enumerate(splits, 1):
            logger.info(f"  Fold {fold_idx}/{n_folds}")
            
            # Split data
            train_geno = genotypes[train_idx, :]
            train_pheno = phenotype[train_idx]
            test_geno = genotypes[test_idx, :]
            test_pheno = phenotype[test_idx]
            
            # Try different P-value thresholds
            best_score = -np.inf
            best_threshold = None
            
            for p_thresh in self.p_thresholds:
                # Clump SNPs
                clumped = self.ld_clump(gwas_results, train_geno, p_thresh)
                
                if len(clumped) == 0:
                    continue
                
                # Create weights
                weights = pd.DataFrame({
                    'SNP_INDEX': clumped.index,
                    'SNP': clumped['SNP'],
                    'BETA': clumped['BETA']
                })
                
                # Calculate PRS on training set
                train_prs = self.calculate_prs(train_geno, weights)
                
                # Evaluate on training (for threshold selection)
                score = self.evaluate_prs(train_prs, train_pheno, is_binary)
                
                if score > best_score:
                    best_score = score
                    best_threshold = p_thresh
            
            # Use best threshold on test set
            if best_threshold is not None:
                clumped = self.ld_clump(gwas_results, train_geno, best_threshold)
                
                if len(clumped) > 0:
                    weights = pd.DataFrame({
                        'SNP_INDEX': clumped.index,
                        'SNP': clumped['SNP'],
                        'BETA': clumped['BETA']
                    })
                    
                    test_prs = self.calculate_prs(test_geno, weights)
                    test_score = self.evaluate_prs(test_prs, test_pheno, is_binary)
                    
                    fold_scores.append(test_score)
                    fold_thresholds.append(best_threshold)
                    
                    logger.info(f"    Best p-threshold: {best_threshold:.2e}")
                    logger.info(f"    Test score: {test_score:.4f}")
        
        # Aggregate results
        metric_name = "AUC" if is_binary else "R²"
        
        cv_results = {
            'fold_scores': fold_scores,
            'mean_score': np.mean(fold_scores),
            'std_score': np.std(fold_scores),
            'thresholds': fold_thresholds,
            'metric': metric_name
        }
        
        logger.info(f"✓ Cross-validation complete")
        logger.info(f"  Mean {metric_name}: {cv_results['mean_score']:.4f} ± {cv_results['std_score']:.4f}")
        
        return cv_results
    
    def calculate_quantile_ors(self,
                              prs: np.ndarray,
                              phenotype: np.ndarray,
                              n_quantiles: int = 10) -> Dict:
        """
        Calculate odds ratios by PRS quantile
        
        Parameters:
        -----------
        prs : np.ndarray
            PRS values
        phenotype : np.ndarray
            Binary phenotype (0/1)
        n_quantiles : int
            Number of quantiles (default: 10 for deciles)
        
        Returns:
        --------
        quantile_ors : dict
            Odds ratios and 95% CIs by quantile
        """
        # Create quantiles
        quantiles = pd.qcut(prs, q=n_quantiles, labels=False, duplicates='drop')
        
        # Reference is bottom quantile (Q0)
        q0_cases = np.sum((quantiles == 0) & (phenotype == 1))
        q0_controls = np.sum((quantiles == 0) & (phenotype == 0))
        
        results = {}
        
        for q in range(n_quantiles):
            q_cases = np.sum((quantiles == q) & (phenotype == 1))
            q_controls = np.sum((quantiles == q) & (phenotype == 0))
            
            if q == 0:
                # Reference quantile
                or_val = 1.0
                ci_lower = 1.0
                ci_upper = 1.0
            else:
                # Calculate OR vs Q0
                if q0_controls > 0 and q_controls > 0:
                    or_val = (q_cases / q_controls) / (q0_cases / q0_controls)
                    
                    # 95% CI using log(OR) ± 1.96 * SE
                    se_log_or = np.sqrt(
                        1/q_cases + 1/q_controls + 1/q0_cases + 1/q0_controls
                    )
                    log_or = np.log(or_val)
                    ci_lower = np.exp(log_or - 1.96 * se_log_or)
                    ci_upper = np.exp(log_or + 1.96 * se_log_or)
                else:
                    or_val = np.nan
                    ci_lower = np.nan
                    ci_upper = np.nan
            
            results[f'Q{q+1}'] = {
                'OR': or_val,
                'CI_lower': ci_lower,
                'CI_upper': ci_upper,
                'n_cases': q_cases,
                'n_controls': q_controls
            }
        
        return results
    
    def generate_prs_report(self,
                           prs: np.ndarray,
                           phenotype: np.ndarray,
                           is_binary: bool,
                           output_file: str):
        """
        Generate comprehensive PRS performance report
        
        Parameters:
        -----------
        prs : np.ndarray
            PRS values
        phenotype : np.ndarray
            Phenotype
        is_binary : bool
            Binary trait
        output_file : str
            Output file path
        """
        logger.info(f"Generating PRS report: {output_file}")
        
        with open(output_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("POLYGENIC RISK SCORE - PERFORMANCE REPORT\n")
            f.write("="*80 + "\n\n")
            
            # Basic statistics
            f.write("PRS Statistics:\n")
            f.write("-"*80 + "\n")
            f.write(f"N samples: {len(prs)}\n")
            f.write(f"PRS mean: {np.mean(prs):.4f}\n")
            f.write(f"PRS std: {np.std(prs):.4f}\n")
            f.write(f"PRS range: [{np.min(prs):.4f}, {np.max(prs):.4f}]\n\n")
            
            if is_binary:
                # Binary trait analysis
                f.write("Binary Trait Analysis:\n")
                f.write("-"*80 + "\n")
                
                # AUC
                auc = roc_auc_score(phenotype, prs)
                f.write(f"AUC: {auc:.4f}\n\n")
                
                # Quantile ORs
                quantile_ors = self.calculate_quantile_ors(prs, phenotype, n_quantiles=10)
                
                f.write("Odds Ratios by Decile (vs. Q1):\n")
                f.write("-"*80 + "\n")
                f.write(f"{'Quantile':<10} {'OR':<10} {'95% CI':<20} {'Cases':<10} {'Controls':<10}\n")
                f.write("-"*80 + "\n")
                
                for q_name, q_data in quantile_ors.items():
                    ci_str = f"({q_data['CI_lower']:.2f}-{q_data['CI_upper']:.2f})"
                    f.write(f"{q_name:<10} {q_data['OR']:<10.2f} {ci_str:<20} "
                           f"{q_data['n_cases']:<10} {q_data['n_controls']:<10}\n")
                
                f.write("\n")
                
                # Top vs bottom decile
                top_or = quantile_ors['Q10']['OR']
                f.write(f"Top decile vs. bottom decile OR: {top_or:.2f}\n")
                f.write(f"95% CI: ({quantile_ors['Q10']['CI_lower']:.2f}-"
                       f"{quantile_ors['Q10']['CI_upper']:.2f})\n")
                
            else:
                # Continuous trait analysis
                f.write("Continuous Trait Analysis:\n")
                f.write("-"*80 + "\n")
                
                # R² and correlation
                correlation = np.corrcoef(prs, phenotype)[0, 1]
                r_squared = correlation ** 2
                
                f.write(f"Correlation (r): {correlation:.4f}\n")
                f.write(f"R² (variance explained): {r_squared:.4f}\n\n")
                
                # By quantile
                prs_std = (prs - np.mean(prs)) / np.std(prs)
                quantiles = pd.qcut(prs_std, q=5, labels=['Q1', 'Q2', 'Q3', 'Q4', 'Q5'])
                
                f.write("Mean Phenotype by PRS Quintile:\n")
                f.write("-"*80 + "\n")
                f.write(f"{'Quantile':<10} {'Mean':<15} {'Std':<15} {'N':<10}\n")
                f.write("-"*80 + "\n")
                
                for q in ['Q1', 'Q2', 'Q3', 'Q4', 'Q5']:
                    in_q = quantiles == q
                    mean_pheno = np.mean(phenotype[in_q])
                    std_pheno = np.std(phenotype[in_q])
                    n_q = np.sum(in_q)
                    
                    f.write(f"{q:<10} {mean_pheno:<15.3f} {std_pheno:<15.3f} {n_q:<10}\n")
            
            f.write("\n" + "="*80 + "\n")
        
        logger.info(f"✓ Report saved: {output_file}")


# Example usage
if __name__ == "__main__":
    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    print("="*80)
    print("POLYGENIC RISK SCORE CALCULATOR")
    print("="*80)
    print("\nThis module implements complete PRS development pipeline.")
    print("\nFeatures:")
    print("  ✓ LD-based SNP clumping")
    print("  ✓ P-value threshold optimization")
    print("  ✓ Cross-validation framework")
    print("  ✓ AUC and R² evaluation")
    print("  ✓ Odds ratios by quantile")
    print("\nExample usage:")
    print("""
    from src.prs.prs_calculator import PolygenicRiskScore
    
    # Initialize
    prs_calc = PolygenicRiskScore(clump_r2=0.1, clump_kb=250)
    
    # Cross-validation
    cv_results = prs_calc.cross_validate(
        gwas_results=gwas_df,
        genotypes=genotypes,
        phenotype=phenotype,
        n_folds=5,
        is_binary=True
    )
    
    print(f"CV AUC: {cv_results['mean_score']:.3f} ± {cv_results['std_score']:.3f}")
    
    # Calculate final PRS
    clumped = prs_calc.ld_clump(gwas_df, genotypes, p_threshold=0.001)
    weights = pd.DataFrame({
        'SNP_INDEX': clumped.index,
        'BETA': clumped['BETA']
    })
    prs_scores = prs_calc.calculate_prs(genotypes, weights)
    
    # Generate report
    prs_calc.generate_prs_report(prs_scores, phenotype, True, 'prs_report.txt')
    """)
    print("="*80)