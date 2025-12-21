"""
Real Data Loader - Pure Python Implementation using scikit-allel

Load and process real genomic data using Python-native packages:
- scikit-allel for VCF files (no external tools needed)
- Works on Windows, Mac, Linux
- No compilation required

Author: Luke M
Repository: https://github.com/lukem2020/gwas-tutorial
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Tuple, Optional, List, Dict
import logging
import urllib.request
import gzip
import yaml

try:
    import allel

    SCIKIT_ALLEL_AVAILABLE = True
except ImportError:
    SCIKIT_ALLEL_AVAILABLE = False
    print("Warning: scikit-allel not installed. Install with: pip install scikit-allel")

logger = logging.getLogger(__name__)

# Known disease-associated SNPs for creating phenotypes
KNOWN_DISEASE_SNPS = {
    "type2diabetes": [
        "rs7903146",  # TCF7L2 - strongest T2D association
        "rs1801282",  # PPARG
        "rs5219",  # KCNJ11
        "rs7754840",  # CDKAL1
        "rs10811661",  # CDKN2A/B
        "rs1111875",  # HHEX
        "rs4402960",  # IGF2BP2
        "rs10010131",  # WFS1
        "rs1801214",  # CAPN10
        "rs13266634",  # SLC30A8
    ],
    "height": [
        "rs6570507",
        "rs17817449",
        "rs6817306",
        "rs1042725",
    ],
    "bmi": [
        "rs9939609",  # FTO
        "rs17782313",  # MC4R
    ],
}


class RealDataLoader:
    """Load and process real genomic data using pure Python packages"""

    def __init__(
        self, config_path: str = "config.yaml", data_dir: Optional[str] = None
    ):
        """
        Initialize data loader with config file

        Parameters:
        -----------
        config_path : str
            Path to config.yaml file
        data_dir : str, optional
            Override data directory from config (for backwards compatibility)
        """
        # Load config file
        config_path = Path(config_path)
        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")

        with open(config_path, "r") as f:
            self.config = yaml.safe_load(f)

        # Get paths from config
        if data_dir is None:
            self.data_dir = Path(self.config["data_dir"]["1000genomes"]).parent
        else:
            self.data_dir = Path(data_dir)

        # Store file paths from config
        self.vcf_path = Path(self.config["file_paths"]["vcf_file"])
        self.panel_path = Path(self.config["file_paths"]["panel_file"])
        self.gwas_sumstats_dir = Path(self.config["data_dir"]["gwas_sumstats"])

        # Create directories if they don't exist
        self.vcf_path.parent.mkdir(parents=True, exist_ok=True)
        self.panel_path.parent.mkdir(parents=True, exist_ok=True)
        self.gwas_sumstats_dir.mkdir(parents=True, exist_ok=True)

        self.logger = logging.getLogger(self.__class__.__name__)

        if not SCIKIT_ALLEL_AVAILABLE:
            raise ImportError(
                "scikit-allel is required but not installed.\n"
                "Install with: pip install scikit-allel"
            )

    def download_1000genomes(self, chromosome: int = 22, force: bool = False):
        """
        Download 1000 Genomes Project data to paths specified in config.yaml

        Parameters:
        -----------
        chromosome : int
            Chromosome number (1-22)
        force : bool
            Re-download even if file exists
        """
        base_url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"

        chr_str = f"chr{chromosome}"
        vcf_filename = f"ALL.{chr_str}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        panel_filename = "integrated_call_samples_v3.20130502.ALL.panel"

        print(f"\n{'='*80}")
        print(f"DOWNLOADING 1000 GENOMES PROJECT DATA - CHROMOSOME {chromosome}")
        print(f"{'='*80}\n")
        print("Note: Using scikit-allel (pure Python, works on Windows)")
        print("No external tools (tabix, bcftools) required!\n")
        print(f"VCF will be saved to: {self.vcf_path}")
        print(f"Panel will be saved to: {self.panel_path}\n")

        # Download VCF to config-specified path
        if self.vcf_path.exists() and not force:
            print(f"✓ VCF file already exists at {self.vcf_path}")
        else:
            print(f"Downloading {vcf_filename}...")
            print(f"Size: ~500MB-2GB (chromosome {chromosome})")
            print(f"URL: {base_url}{vcf_filename}")

            try:
                urllib.request.urlretrieve(
                    f"{base_url}{vcf_filename}",
                    self.vcf_path,
                    reporthook=self._download_progress,
                )
                print(f"\n✓ Downloaded VCF to {self.vcf_path}")
            except Exception as e:
                raise RuntimeError(f"Download failed: {e}")

        # Download panel to config-specified path
        if self.panel_path.exists() and not force:
            print(f"✓ Panel file already exists at {self.panel_path}")
        else:
            print(f"\nDownloading {panel_filename}...")
            urllib.request.urlretrieve(f"{base_url}{panel_filename}", self.panel_path)
            print(f"✓ Downloaded panel to {self.panel_path}")

        print(f"\n✓ 1000 Genomes chromosome {chromosome} data ready!")
        print(f"VCF location: {self.vcf_path}")
        print(f"Panel location: {self.panel_path}")

        return self.vcf_path.parent

    def _download_progress(self, block_num, block_size, total_size):
        """Progress callback for downloads"""
        downloaded = block_num * block_size
        if total_size > 0:
            percent = min(downloaded * 100.0 / total_size, 100)
            mb_downloaded = downloaded / (1024 * 1024)
            mb_total = total_size / (1024 * 1024)

            # Print progress every 5%
            if block_num % 100 == 0:
                print(
                    f"\rProgress: {percent:.1f}% ({mb_downloaded:.1f}/{mb_total:.1f} MB)",
                    end="",
                )

    def load_1000genomes(
        self,
        chromosome: int = 22,
        max_snps: Optional[int] = None,
        populations: List[str] = ["EUR"],
        maf: float = 0.01,
        region: Optional[str] = None,
    ) -> Tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
        """
        Load 1000 Genomes VCF using scikit-allel (pure Python)

        Parameters:
        -----------
        chromosome : int
            Chromosome to load
        max_snps : int, optional
            Limit SNPs (for testing)
        populations : list
            Population codes: EUR, AFR, EAS, SAS, AMR, ALL
        maf : float
            Minimum minor allele frequency
        region : str, optional
            Genomic region (e.g., "22:16000000-17000000")

        Returns:
        --------
        genotypes : np.ndarray (n_samples x n_snps)
            Genotype matrix (0, 1, 2 = number of alt alleles)
        snp_info : pd.DataFrame (CHR, SNP, BP, REF, ALT)
        sample_info : pd.DataFrame (sample ID, population info)
        """
        self.logger.info(f"Loading 1000 Genomes chr{chromosome}...")

        # Use paths from config
        vcf_file = self.vcf_path
        panel_file = self.panel_path

        if not vcf_file.exists():
            raise FileNotFoundError(
                f"1000 Genomes VCF file not found at {vcf_file}.\n"
                f"Run: loader.download_1000genomes(chromosome={chromosome})"
            )

        if not panel_file.exists():
            raise FileNotFoundError(
                f"1000 Genomes panel file not found at {panel_file}.\n"
                f"Run: loader.download_1000genomes(chromosome={chromosome})"
            )

        # Load panel (population information)
        panel = pd.read_csv(panel_file, sep="\t")

        # Filter by population
        if "ALL" not in populations:
            panel = panel[panel["super_pop"].isin(populations)]

        samples_to_keep = panel["sample"].tolist()

        print(f"\nLoading VCF using scikit-allel (pure Python)...")
        print(f"Chromosome: {chromosome}")
        print(f"Populations: {populations}")
        print(f"Samples: {len(samples_to_keep)}")
        if region:
            print(f"Region: {region}")

        # Read VCF file
        print("Reading VCF file (this may take a few minutes)...")

        fields = [
            "variants/CHROM",
            "variants/POS",
            "variants/ID",
            "variants/REF",
            "variants/ALT",
            "calldata/GT",
            "samples",
        ]

        if region:
            callset = allel.read_vcf(str(vcf_file), fields=fields, region=region)
        else:
            callset = allel.read_vcf(str(vcf_file), fields=fields)

        if callset is None:
            raise ValueError("Failed to read VCF file")

        print(
            f"✓ Loaded VCF: {len(callset['variants/POS'])} variants, {len(callset['samples'])} samples"
        )

        # Get sample indices for selected populations
        all_samples = callset["samples"]
        sample_indices = [i for i, s in enumerate(all_samples) if s in samples_to_keep]

        # Extract genotypes for selected samples
        gt = callset["calldata/GT"][:, sample_indices, :]

        # Convert to allele counts (0, 1, 2)
        # Sum across two chromosomes
        genotypes = gt.sum(axis=2)

        # Handle missing data (-1 becomes NaN)
        genotypes = genotypes.astype(float)
        genotypes[genotypes < 0] = np.nan

        # Transpose to (samples x SNPs)
        genotypes = genotypes.T

        # Create SNP info dataframe
        snp_info = pd.DataFrame(
            {
                "CHR": callset["variants/CHROM"],
                "SNP": callset["variants/ID"],
                "BP": callset["variants/POS"],
                "REF": [
                    ref.decode() if isinstance(ref, bytes) else str(ref)
                    for ref in callset["variants/REF"]
                ],
                "ALT": [
                    alt[0].decode() if isinstance(alt[0], bytes) else str(alt[0])
                    for alt in callset["variants/ALT"]
                ],
            }
        )

        # Filter by MAF
        if maf > 0:
            print(f"\nFiltering by MAF > {maf}...")
            af = np.nanmean(genotypes, axis=0) / 2
            maf_values = np.minimum(af, 1 - af)
            maf_pass = maf_values >= maf

            genotypes = genotypes[:, maf_pass]
            snp_info = snp_info[maf_pass].reset_index(drop=True)

            print(f"✓ Retained {genotypes.shape[1]} SNPs after MAF filter")

        # Limit SNPs if requested
        if max_snps and max_snps < genotypes.shape[1]:
            print(f"\nLimiting to {max_snps} SNPs for testing...")
            genotypes = genotypes[:, :max_snps]
            snp_info = snp_info.iloc[:max_snps].reset_index(drop=True)

        # Create sample info
        sample_info = pd.DataFrame({"IID": [all_samples[i] for i in sample_indices]})
        sample_info = sample_info.merge(
            panel[["sample", "pop", "super_pop"]],
            left_on="IID",
            right_on="sample",
            how="left",
        )

        print(
            f"\n✓ Final dataset: {genotypes.shape[0]} samples x {genotypes.shape[1]} SNPs"
        )

        return genotypes, snp_info, sample_info

    def load_gwas_sumstats(
        self, file_path: str, sep: str = r"\s+", compression: str = "infer"
    ) -> pd.DataFrame:
        """
        Load GWAS summary statistics from public sources

        Parameters:
        -----------
        file_path : str
            Path to summary statistics file (relative paths will be resolved
            against config's gwas_sumstats directory)
        sep : str
            Column separator
        compression : str
            'infer', 'gzip', None

        Returns:
        --------
        sumstats : pd.DataFrame
            Standardized columns: SNP, CHR, BP, A1, A2, BETA, SE, P
        """
        self.logger.info(f"Loading GWAS summary statistics: {file_path}")

        # Resolve file path - if relative, try in gwas_sumstats directory from config
        file_path_obj = Path(file_path)
        if not file_path_obj.is_absolute() and not file_path_obj.exists():
            # Try in gwas_sumstats directory from config
            alt_path = self.gwas_sumstats_dir / file_path
            if alt_path.exists():
                file_path = str(alt_path)
                print(f"Found file in config directory: {file_path}")

        # Read file
        df = pd.read_csv(file_path, sep=sep, compression=compression, low_memory=False)

        print(f"Original columns: {df.columns.tolist()}")

        # Standardize column names
        column_mapping = {
            # SNP ID
            "rsid": "SNP",
            "rs": "SNP",
            "snp": "SNP",
            "MarkerName": "SNP",
            "variant_id": "SNP",
            "ID": "SNP",
            # Chromosome
            "chr": "CHR",
            "chromosome": "CHR",
            "#chr": "CHR",
            # Position
            "pos": "BP",
            "position": "BP",
            "bp": "BP",
            "base_pair_location": "BP",
            # Effect size
            "effect": "BETA",
            "beta": "BETA",
            "b": "BETA",
            "Effect": "BETA",
            # Standard error
            "se": "SE",
            "stderr": "SE",
            "standard_error": "SE",
            "StdErr": "SE",
            # P-value
            "p": "P",
            "pval": "P",
            "p_value": "P",
            "pvalue": "P",
            "P-value": "P",
            "p.value": "P",
            "P_VALUE": "P",
            # Alleles
            "A1": "A1",
            "effect_allele": "A1",
            "Allele1": "A1",
            "EA": "A1",
            "A2": "A2",
            "other_allele": "A2",
            "Allele2": "A2",
            "NEA": "A2",
            "ref": "A2",
            "reference_allele": "A2",
            # Frequency
            "eaf": "EAF",
            "freq": "EAF",
            "effect_allele_freq": "EAF",
            "Freq1": "EAF",
            # Odds ratio
            "OR": "OR",
            "odds_ratio": "OR",
            "or": "OR",
            # Sample size
            "n": "N",
            "N_total": "N",
            "TotalSampleSize": "N",
        }

        # Apply mapping
        df = df.rename(
            columns={k: v for k, v in column_mapping.items() if k in df.columns}
        )

        # Convert OR to BETA if needed
        if "OR" in df.columns and "BETA" not in df.columns:
            df["BETA"] = np.log(df["OR"])
            print("✓ Converted OR to BETA")

        # Calculate SE from P and BETA if missing
        if "SE" not in df.columns and "BETA" in df.columns and "P" in df.columns:
            from scipy import stats

            z = stats.norm.ppf(1 - df["P"] / 2)
            z = np.where(z == 0, np.nan, z)
            df["SE"] = np.abs(df["BETA"] / z)
            print("✓ Calculated SE from P and BETA")

        # Remove rows with missing critical values
        before = len(df)
        df = df.dropna(subset=["SNP", "P"])
        after = len(df)
        if before != after:
            print(f"Removed {before - after} SNPs with missing SNP ID or P-value")

        # Remove invalid P-values
        df = df[(df["P"] > 0) & (df["P"] <= 1)]

        print(f"\n✓ Loaded {len(df):,} SNPs from summary statistics")
        print(f"Final columns: {df.columns.tolist()}")

        return df

    def create_simulated_phenotype(
        self,
        genotypes: np.ndarray,
        snp_info: pd.DataFrame,
        trait: str = "type2diabetes",
        heritability: float = 0.3,
        binary: bool = True,
        case_prevalence: float = 0.5,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Create phenotype using known disease SNPs

        Parameters:
        -----------
        genotypes : np.ndarray
            Genotype matrix
        snp_info : pd.DataFrame
            SNP information
        trait : str
            Disease/trait name
        heritability : float
            Proportion of variance explained
        binary : bool
            Create binary phenotype
        case_prevalence : float
            Proportion of cases

        Returns:
        --------
        phenotype : np.ndarray
            Simulated phenotype
        genetic_value : np.ndarray
            True genetic component
        """
        if trait not in KNOWN_DISEASE_SNPS:
            raise ValueError(
                f"Unknown trait. Available: {list(KNOWN_DISEASE_SNPS.keys())}"
            )

        known_snps = KNOWN_DISEASE_SNPS[trait]

        # Find which SNPs are in our data
        snp_indices = []
        found_snps = []

        for snp in known_snps:
            if snp in snp_info["SNP"].values:
                idx = snp_info[snp_info["SNP"] == snp].index[0]
                snp_indices.append(idx)
                found_snps.append(snp)

        if len(snp_indices) == 0:
            raise ValueError(f"No known {trait} SNPs found in genotype data")

        print(f"\nCreating {trait} phenotype:")
        print(f"  Using {len(found_snps)}/{len(known_snps)} known {trait} SNPs")
        print(f"  SNPs found: {found_snps}")

        # Calculate genetic value
        np.random.seed(42)
        effects = np.random.normal(0, 0.3, len(snp_indices))

        genetic_value = genotypes[:, snp_indices] @ effects

        # Add environmental noise
        var_genetic = np.var(genetic_value)
        var_environmental = var_genetic * (1 - heritability) / heritability
        environmental_noise = np.random.normal(
            0, np.sqrt(var_environmental), len(genetic_value)
        )

        phenotype_continuous = genetic_value + environmental_noise

        if binary:
            threshold = np.quantile(phenotype_continuous, 1 - case_prevalence)
            phenotype = (phenotype_continuous > threshold).astype(int)

            n_cases = np.sum(phenotype == 1)
            n_controls = np.sum(phenotype == 0)

            print(f"  Cases: {n_cases}")
            print(f"  Controls: {n_controls}")
            print(f"  Prevalence: {n_cases/(n_cases+n_controls):.2%}")
        else:
            phenotype = phenotype_continuous
            print(f"  Mean: {phenotype.mean():.3f}")
            print(f"  Std: {phenotype.std():.3f}")

        return phenotype, genetic_value


def quick_load_real_data(
    chromosome: int = 22,
    populations: List[str] = ["EUR"],
    max_snps: Optional[int] = 10000,
    trait: str = "type2diabetes",
    region: Optional[str] = None,
) -> Tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """
    Quick function to load real data for testing

    Example:
    --------
    genotypes, phenotype, snp_info = quick_load_real_data(
        chromosome=22,
        max_snps=10000,
        trait='type2diabetes'
    )
    """
    loader = RealDataLoader()

    print("Loading 1000 Genomes data...")
    genotypes, snp_info, sample_info = loader.load_1000genomes(
        chromosome=chromosome, populations=populations, max_snps=max_snps, region=region
    )

    print("\nCreating simulated phenotype from known disease SNPs...")
    phenotype, genetic_value = loader.create_simulated_phenotype(
        genotypes, snp_info, trait=trait
    )

    return genotypes, phenotype, snp_info


if __name__ == "__main__":
    import sys

    print("Real Data Loader - Pure Python Implementation")
    print("=" * 80)
    print("\nFeatures:")
    print("✓ Uses scikit-allel (pure Python)")
    print("✓ Works on Windows, Mac, Linux")
    print("✓ No external tools required")
    print("✓ No compilation needed")

    # Initialize loader
    try:
        loader = RealDataLoader()
        print(f"\n✓ Config loaded from config.yaml")
        print(f"  VCF path: {loader.vcf_path}")
        print(f"  Panel path: {loader.panel_path}")
        print(f"  GWAS sumstats dir: {loader.gwas_sumstats_dir}")
    except Exception as e:
        print(f"\n✗ Error loading config: {e}")
        sys.exit(1)

    # Check if data exists
    vcf_exists = loader.vcf_path.exists()
    panel_exists = loader.panel_path.exists()

    print(f"\n{'='*80}")
    print("Data Status:")
    print(
        f"  VCF file: {'✓ EXISTS' if vcf_exists else '✗ MISSING'} - {loader.vcf_path}"
    )
    print(
        f"  Panel file: {'✓ EXISTS' if panel_exists else '✗ MISSING'} - {loader.panel_path}"
    )

    # Automatically download missing files
    if not vcf_exists or not panel_exists:
        print(f"\n{'='*80}")
        print("DOWNLOADING MISSING FILES...")
        print(f"{'='*80}")
        try:
            # Extract chromosome from VCF filename if possible, default to 22
            chromosome = 22
            vcf_name = loader.vcf_path.name
            if "chr" in vcf_name.lower():
                # Try to extract chromosome number from filename
                import re

                match = re.search(r"chr(\d+)", vcf_name, re.IGNORECASE)
                if match:
                    chromosome = int(match.group(1))

            loader.download_1000genomes(chromosome=chromosome, force=False)
            print("\n✓ Download complete!")

            # Verify files now exist
            vcf_exists = loader.vcf_path.exists()
            panel_exists = loader.panel_path.exists()
        except Exception as e:
            print(f"\n✗ Download failed: {e}")
            import traceback

            traceback.print_exc()
            sys.exit(1)

    if vcf_exists and panel_exists:
        print(f"\n{'='*80}")
        print("✓ All required data files are available!")
        print("\nTo load data, use:")
        print(
            "  from src.utils.read_data_loader import RealDataLoader, quick_load_real_data"
        )
        print("  loader = RealDataLoader()")
        print(
            "  genotypes, snp_info, sample_info = loader.load_1000genomes(chromosome=22, max_snps=1000)"
        )
        print("\nOr use the quick function:")
        print("  genotypes, phenotype, snp_info = quick_load_real_data(max_snps=1000)")

    # Handle command line arguments
    if len(sys.argv) > 1:
        if sys.argv[1] == "download":
            print(f"\n{'='*80}")
            print("DOWNLOADING 1000 GENOMES DATA")
            print(f"{'='*80}")
            try:
                loader.download_1000genomes(chromosome=22, force=False)
                print("\n✓ Download complete!")
            except Exception as e:
                print(f"\n✗ Download failed: {e}")
                sys.exit(1)
        elif sys.argv[1] == "demo":
            print(f"\n{'='*80}")
            print("RUNNING DEMO")
            print(f"{'='*80}")
            # Ensure files exist before demo
            if not loader.vcf_path.exists() or not loader.panel_path.exists():
                print("Downloading missing files first...")
                chromosome = 22
                vcf_name = loader.vcf_path.name
                if "chr" in vcf_name.lower():
                    import re

                    match = re.search(r"chr(\d+)", vcf_name, re.IGNORECASE)
                    if match:
                        chromosome = int(match.group(1))
                loader.download_1000genomes(chromosome=chromosome, force=False)
            try:
                print("\nLoading small subset of data (1000 SNPs)...")
                genotypes, phenotype, snp_info = quick_load_real_data(
                    chromosome=22, max_snps=1000, populations=["EUR"]
                )
                print(f"\n✓ Demo complete!")
                print(f"  Genotypes shape: {genotypes.shape}")
                print(f"  Phenotype shape: {phenotype.shape}")
                print(f"  SNPs loaded: {len(snp_info)}")
            except Exception as e:
                print(f"\n✗ Demo failed: {e}")
                import traceback

                traceback.print_exc()
                sys.exit(1)
        else:
            print(f"\nUnknown command: {sys.argv[1]}")
            print("Available commands: download, demo")
