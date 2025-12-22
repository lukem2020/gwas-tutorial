"""
Download all GEO datasets listed in config.yaml.
"""

from src.download_geo_data import download_all_geo_datasets
import sys


def main():
    """Download all GEO datasets from config."""
    print("=" * 60)
    print("Download GEO Datasets")
    print("=" * 60)
    print()
    
    try:
        results = download_all_geo_datasets()
        
        print()
        print("=" * 60)
        print("Download Summary")
        print("=" * 60)
        successful = [gse for gse, data in results.items() if 'error' not in data]
        failed = [gse for gse, data in results.items() if 'error' in data]
        
        print(f"Successfully downloaded: {len(successful)}/{len(results)}")
        if successful:
            for gse in successful:
                data = results[gse]
                expr_shape = data['expression'].shape
                n_samples = len(data['phenotypes'])
                print(f"  ✓ {gse}: {expr_shape[0]} genes × {expr_shape[1]} samples, {n_samples} phenotypes")
        
        if failed:
            print(f"\nFailed: {len(failed)}")
            for gse in failed:
                print(f"  ✗ {gse}: {results[gse].get('error', 'Unknown error')}")
        
        print()
        print("Files saved to: data/expression/")
        print("  - {GSE}_expression.csv (expression data)")
        print("  - {GSE}_phenotypes.csv (phenotype data)")
        
        return 0 if not failed else 1
        
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())

