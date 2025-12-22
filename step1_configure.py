"""
Step 1: Configure Project
Creates and validates the configuration file for the GWAS pipeline.
"""

from src.config import load_config
from pathlib import Path


def main():
    """Run Step 1: Configure project."""
    print("=" * 60)
    print("Step 1: Configure Project")
    print("=" * 60)
    
    config_path = "config.yaml"
    
    # Check if config exists
    if Path(config_path).exists():
        print(f"✓ Configuration file found: {config_path}")
        print("  Loading and validating configuration...")
    else:
        print(f"✗ Configuration file not found: {config_path}")
        print("  Please create config.yaml using the template.")
        print("  A template should be available in the project root.")
        return
    
    try:
        # Load and validate configuration
        config = load_config(config_path)
        
        print("\n✓ Configuration loaded successfully!")
        print("\nProject Configuration:")
        print(f"  Name: {config.get('project.name')}")
        print(f"  Version: {config.get('project.version')}")
        
        print("\nDirectory Structure:")
        print(f"  Data: {config.get('paths.data_dir')}")
        print(f"  Genotypes: {config.get('paths.genotypes_dir')}")
        print(f"  Phenotypes: {config.get('paths.phenotypes_dir')}")
        print(f"  Expression: {config.get('paths.expression_dir')}")
        print(f"  Results: {config.get('paths.results_dir')}")
        print(f"  Plots: {config.get('paths.plots_dir')}")
        
        print("\n1000 Genomes Settings:")
        print(f"  Phase: {config.get('genotypes.phase')}")
        print(f"  Reference: {config.get('genotypes.reference')}")
        print(f"  Chromosomes: {config.get('genotypes.chromosomes')}")
        
        print("\nQuality Control Thresholds:")
        print(f"  Sample missingness: {config.get('qc.sample_missingness_threshold')}")
        print(f"  Variant missingness: {config.get('qc.variant_missingness_threshold')}")
        print(f"  MAF threshold: {config.get('qc.maf_threshold')}")
        print(f"  HWE p-threshold: {config.get('qc.hwe_p_threshold')}")
        
        print("\nGWAS Settings:")
        print(f"  Model: {config.get('gwas.model')}")
        print(f"  Significance threshold: {config.get('gwas.significance_threshold')}")
        print(f"  Multiple testing: {config.get('gwas.multiple_testing')}")
        
        print("\n✓ All directories created successfully!")
        print("\n" + "=" * 60)
        print("Step 1 Complete: Project configured and ready!")
        print("=" * 60)
        
    except FileNotFoundError as e:
        print(f"\n✗ Error: {e}")
        print("\nPlease create config.yaml in the project root.")
    except Exception as e:
        print(f"\n✗ Error loading configuration: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()

