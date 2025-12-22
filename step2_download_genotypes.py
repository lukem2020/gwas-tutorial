"""
Step 2: Download Genotype Data (1000 Genomes)
Downloads VCF files for specified chromosomes from 1000 Genomes Project.
"""

from src.config import load_config
from src.download_genotypes import download_genotypes, GenotypeDownloader
from pathlib import Path
import sys


def main():
    """Run Step 2: Download genotype data."""
    print("=" * 60)
    print("Step 2: Download Genotype Data (1000 Genomes)")
    print("=" * 60)
    print()
    
    # Load configuration
    try:
        config = load_config()
        print("✓ Configuration loaded")
    except Exception as e:
        print(f"✗ Error loading configuration: {e}")
        print("  Please run Step 1 first: python step1_configure.py")
        return 1
    
    # Get settings from config
    output_dir = config.get('paths.genotypes_dir')
    chromosomes = config.get('genotypes.chromosomes')
    phase = config.get('genotypes.phase')
    reference = config.get('genotypes.reference')
    
    print(f"Settings:")
    print(f"  Phase: {phase}")
    print(f"  Reference: {reference}")
    print(f"  Chromosomes: {chromosomes}")
    print(f"  Output directory: {output_dir}")
    print()
    
    # Check if output directory exists
    output_path = Path(output_dir)
    if not output_path.exists():
        print(f"✗ Output directory does not exist: {output_dir}")
        print("  Please run Step 1 first to create directories")
        return 1
    
    # Ask for confirmation if files already exist
    downloader = GenotypeDownloader(output_dir, chromosomes)
    existing = [chr for chr in chromosomes if downloader.is_downloaded(chr)]
    
    if existing:
        print(f"Note: {len(existing)} chromosome file(s) already downloaded:")
        for chr in existing:
            file_size_mb = downloader.get_local_path(chr).stat().st_size / (1024 * 1024)
            print(f"  - Chromosome {chr} ({file_size_mb:.1f} MB)")
        print()
        response = input("Re-download existing files? (y/N): ").strip().lower()
        force = response == 'y'
    else:
        force = False
    
    print()
    print("Starting download...")
    print("(This may take 10-60 minutes depending on your internet connection)")
    print()
    
    # Download genotypes
    try:
        results = download_genotypes(config, force=force)
        
        # Print summary
        print()
        print("=" * 60)
        print("Download Summary")
        print("=" * 60)
        print(f"  Successfully downloaded: {len(results['success'])} chromosome(s)")
        if results['success']:
            print(f"    Chromosomes: {results['success']}")
        
        print(f"  Already existed (skipped): {len(results['skipped'])} chromosome(s)")
        if results['skipped']:
            print(f"    Chromosomes: {results['skipped']}")
        
        if results['failed']:
            print(f"  Failed: {len(results['failed'])} chromosome(s)")
            print(f"    Chromosomes: {results['failed']}")
            print()
            print("  You may need to:")
            print("    - Check your internet connection")
            print("    - Verify the 1000 Genomes FTP server is accessible")
            print("    - Try running again later")
            return 1
        
        # Verify downloads
        print()
        print("Verifying downloads...")
        downloader = GenotypeDownloader(output_dir, chromosomes)
        verification = downloader.verify_downloads()
        
        if len(verification['valid']) == len(chromosomes):
            print(f"✓ All {len(verification['valid'])} chromosome files verified")
        else:
            print(f"⚠ Warning: {len(verification['invalid'])} file(s) missing or invalid")
            if verification['invalid']:
                print(f"  Invalid: {verification['invalid']}")
        
        print()
        print("=" * 60)
        print("Step 2 Complete: Genotype data downloaded!")
        print("=" * 60)
        print()
        print("Next step: Step 3 - Prepare phenotype data")
        print("  Run: python step3_prepare_phenotypes.py")
        
        return 0
        
    except KeyboardInterrupt:
        print()
        print("\n✗ Download interrupted by user")
        return 1
    except Exception as e:
        print(f"\n✗ Error during download: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())

