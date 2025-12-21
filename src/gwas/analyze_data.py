"""
Quick Start: Analyze Real 1000 Genomes Data

This script demonstrates the complete GWAS pipeline with real data.
Uses only Python packages - works on Windows, Mac, Linux!

Author: Luke M
Repository: https://github.com/lukem2020/gwas-tutorial
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent))

print("=" * 80)
print("GWAS TUTORIAL - ANALYZE REAL 1000 GENOMES DATA")
print("=" * 80)
print("\nThis script will:")
print("  1. Load real genotypes from 1000 Genomes Project")
print("  2. Create phenotype using known disease SNPs")
print("  3. Run GWAS with Linear Mixed Models")
print("  4. Calculate Polygenic Risk Scores")
print("  5. Generate visualizations")
print("\n" + "=" * 80 + "\n")

# =============================================================================
# STEP 1: Load Real Data
# =============================================================================

print("[Step 1/5] Loading Real 1000 Genomes Data...")
print("  This uses scikit-allel (pure Python, works on Windows)")

try:
    from src.utils.read_data_loader import RealDataLoader, quick_load_real_data

    # Load data (start with small subset for quick testing)
    genotypes, phenotype, snp_info = quick_load_real_data(
        chromosome=22,
        populations=["EUR"],
        max_snps=5000,  # Start small - increase to 50000 for full analysis
        trait="type2diabetes",
    )

    print(f"\n✓ Data loaded successfully!")
    print(f"  Samples: {genotypes.shape[0]}")
    print(f"  SNPs: {genotypes.shape[1]}")
    print(f"  Cases: {sum(phenotype == 1)}")
    print(f"  Controls: {sum(phenotype == 0)}")

except ImportError as e:
    print(f"\n✗ Error: {e}")
    print("\nMake sure you have:")
    print("  1. Installed all packages: pip install -r requirements.txt")
    print("  2. Added real_data_loader.py to src/utils/")
    sys.exit(1)

except FileNotFoundError:
    print("\n✗ 1000 Genomes data not found!")
    print("\nDownload it first:")
    print(
        '  python -c "from src.utils.real_data_loader import RealDataLoader; RealDataLoader().download_1000genomes(22)"'
    )
    sys.exit(1)

# =============================================================================
# STEP 2: Basic Association Test (if LMM not available)
# =============================================================================

print("\n[Step 2/5] Running Association Analysis...")

try:
    from src.gwas.linear_mixed_model import LinearMixedModelGWAS

    print("  Using Linear Mixed Model (advanced)")

    lmm = LinearMixedModelGWAS(use_reml=True, min_maf=0.01)
    results = lmm.run_gwas(genotypes=genotypes, phenotype=phenotype, snp_info=snp_info)

    print(f"\n✓ GWAS complete!")
    print(f"  Heritability (h²): {results.heritability:.4f}")
    print(f"  Genomic inflation (λ): {results.genomic_control:.4f}")

    # Get results
    pvalues = results.pvalues
    beta = results.beta

except ImportError:
    print("  Linear Mixed Model not available, using simple logistic regression")

    from scipy import stats
    from sklearn.linear_model import LogisticRegression

    # Simple SNP-by-SNP test
    pvalues = []
    beta_values = []

    print(f"  Testing {genotypes.shape[1]} SNPs...")

    for i in range(genotypes.shape[1]):
        if i % 1000 == 0:
            print(f"    Progress: {i}/{genotypes.shape[1]}")

        snp = genotypes[:, i]

        # Skip if too many missing
        if np.sum(np.isnan(snp)) > len(snp) * 0.1:
            pvalues.append(np.nan)
            beta_values.append(np.nan)
            continue

        # Impute missing with mean
        snp_clean = snp.copy()
        snp_clean[np.isnan(snp)] = np.nanmean(snp)

        # Logistic regression
        try:
            lr = LogisticRegression(penalty=None, max_iter=100)
            lr.fit(snp_clean.reshape(-1, 1), phenotype)

            beta_values.append(lr.coef_[0][0])

            # Calculate p-value using Wald test
            from sklearn.metrics import log_loss

            null_pred = np.full(len(phenotype), phenotype.mean())
            alt_pred = lr.predict_proba(snp_clean.reshape(-1, 1))[:, 1]

            # Simple chi-square test
            chi2, p = stats.chi2_contingency(
                [
                    [
                        sum((phenotype == 0) & (snp_clean == 0)),
                        sum((phenotype == 0) & (snp_clean > 0)),
                    ],
                    [
                        sum((phenotype == 1) & (snp_clean == 0)),
                        sum((phenotype == 1) & (snp_clean > 0)),
                    ],
                ]
            )[:2]

            pvalues.append(p)

        except:
            pvalues.append(np.nan)
            beta_values.append(np.nan)

    pvalues = np.array(pvalues)
    beta = np.array(beta_values)

    print("\n✓ Association testing complete!")

# =============================================================================
# STEP 3: Find Top Hits
# =============================================================================

print("\n[Step 3/5] Analyzing Results...")

# Top associations
valid_idx = ~np.isnan(pvalues)
if np.sum(valid_idx) > 0:

    top_5_idx = np.argsort(pvalues[valid_idx])[:5]
    actual_indices = np.where(valid_idx)[0][top_5_idx]

    print("\nTop 5 Associations:")
    print(
        f"{'Rank':<6} {'SNP':<15} {'Chr':<5} {'Position':<12} {'P-value':<15} {'Beta':<10}"
    )
    print("-" * 70)

    for rank, idx in enumerate(actual_indices, 1):
        print(
            f"{rank:<6} {snp_info.iloc[idx]['SNP']:<15} "
            f"{snp_info.iloc[idx]['CHR']:<5} "
            f"{snp_info.iloc[idx]['BP']:<12} "
            f"{pvalues[idx]:<15.2e} "
            f"{beta[idx]:<10.4f}"
        )

    # Count significant
    n_sig = np.sum(pvalues < 5e-8)
    n_sug = np.sum(pvalues < 1e-5)

    print(f"\nSignificance Summary:")
    print(f"  Genome-wide (p < 5×10⁻⁸): {n_sig}")
    print(f"  Suggestive (p < 1×10⁻⁵): {n_sug}")

# =============================================================================
# STEP 4: Calculate PRS (Simple Version)
# =============================================================================

print("\n[Step 4/5] Calculating Polygenic Risk Score...")

# Simple PRS using p < 0.01 SNPs
sig_snps = pvalues < 0.01
sig_indices = np.where(sig_snps)[0]

if len(sig_indices) > 0:
    print(f"  Using {len(sig_indices)} SNPs with p < 0.01")

    # Calculate PRS
    prs_scores = genotypes[:, sig_indices] @ beta[sig_indices]

    # Standardize
    prs_std = (prs_scores - np.mean(prs_scores)) / np.std(prs_scores)

    # Calculate AUC
    from sklearn.metrics import roc_auc_score

    auc = roc_auc_score(phenotype, prs_std)

    print(f"  PRS AUC: {auc:.4f}")

    # Odds ratios by quantile
    quantiles = pd.qcut(prs_std, q=5, labels=["Q1", "Q2", "Q3", "Q4", "Q5"])

    print("\n  Odds Ratios by PRS Quintile (vs Q1):")
    for q in ["Q1", "Q2", "Q3", "Q4", "Q5"]:
        in_q = quantiles == q
        cases_in_q = np.sum(phenotype[in_q] == 1)
        controls_in_q = np.sum(phenotype[in_q] == 0)

        if q == "Q1":
            cases_ref = cases_in_q
            controls_ref = controls_in_q
            or_val = 1.0
        else:
            or_val = (cases_in_q / controls_in_q) / (cases_ref / controls_ref)

        print(
            f"    {q}: OR = {or_val:.2f} ({cases_in_q} cases, {controls_in_q} controls)"
        )
else:
    print("  No significant SNPs for PRS")

# =============================================================================
# STEP 5: Save Results
# =============================================================================

print("\n[Step 5/5] Saving Results...")

# Create results directory
Path("results").mkdir(exist_ok=True)

# Save GWAS results
gwas_df = pd.DataFrame(
    {
        "CHR": snp_info["CHR"],
        "SNP": snp_info["SNP"],
        "BP": snp_info["BP"],
        "BETA": beta,
        "P": pvalues,
    }
)

gwas_df.to_csv("results/gwas_results.csv", index=False)
print("  ✓ Saved: results/gwas_results.csv")

# Save PRS if calculated
if len(sig_indices) > 0:
    prs_df = pd.DataFrame(
        {
            "IID": range(len(prs_scores)),
            "PRS": prs_scores,
            "PRS_STD": prs_std,
            "PHENOTYPE": phenotype,
        }
    )
    prs_df.to_csv("results/prs_scores.csv", index=False)
    print("  ✓ Saved: results/prs_scores.csv")

# Create simple plot
try:
    import matplotlib.pyplot as plt

    # Manhattan plot (simplified)
    plt.figure(figsize=(14, 6))
    plt.scatter(range(len(pvalues)), -np.log10(pvalues), s=5, alpha=0.6)
    plt.axhline(y=-np.log10(5e-8), color="r", linestyle="--", label="Genome-wide")
    plt.axhline(y=-np.log10(1e-5), color="b", linestyle="--", label="Suggestive")
    plt.xlabel("SNP Index")
    plt.ylabel("-log₁₀(P-value)")
    plt.title("Manhattan Plot (Chromosome 22)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("results/manhattan_plot.png", dpi=300)
    print("  ✓ Saved: results/manhattan_plot.png")
    plt.close()

    # QQ plot
    from scipy import stats

    pvals_clean = pvalues[~np.isnan(pvalues)]
    pvals_sorted = np.sort(pvals_clean)
    expected = np.arange(1, len(pvals_sorted) + 1) / (len(pvals_sorted) + 1)

    plt.figure(figsize=(8, 8))
    plt.scatter(-np.log10(expected), -np.log10(pvals_sorted), alpha=0.5, s=10)
    plt.plot(
        [0, -np.log10(expected).max()],
        [0, -np.log10(expected).max()],
        "r--",
        label="y=x",
    )
    plt.xlabel("Expected -log₁₀(P)")
    plt.ylabel("Observed -log₁₀(P)")
    plt.title("Q-Q Plot")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("results/qq_plot.png", dpi=300)
    print("  ✓ Saved: results/qq_plot.png")
    plt.close()

except Exception as e:
    print(f"  ⚠ Could not create plots: {e}")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

print("\n" + "=" * 80)
print("✓ ANALYSIS COMPLETE!")
print("=" * 80)
print("\nResults saved to results/ directory:")
print("  - gwas_results.csv: All SNP associations")
print("  - prs_scores.csv: Individual risk scores")
print("  - manhattan_plot.png: Genome-wide visualization")
print("  - qq_plot.png: QC plot")
print("\nThis analysis used REAL genomic data from 1000 Genomes Project!")
print("\nNext steps:")
print("  1. Increase max_snps for more comprehensive analysis")
print("  2. Try different populations")
print("  3. Implement full Linear Mixed Models")
print("  4. Add more sophisticated PRS methods")
print("\nRepository: https://github.com/lukem2020/gwas-tutorial")
print("=" * 80)
