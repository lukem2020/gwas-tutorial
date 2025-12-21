"""Utility script to work with VCF files using scikit-allel (Python-only, no external tools needed)."""
import sys
from pathlib import Path
import logging

try:
    import allel
    SCIKIT_ALLEL_AVAILABLE = True
except ImportError:
    SCIKIT_ALLEL_AVAILABLE = False

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())


def read_vcf(vcf_path: str, region: str = None):
    """
    Read a VCF file using scikit-allel (Python-only, works on Windows).
    
    Args:
        vcf_path: Path to the VCF file (.vcf or .vcf.gz)
        region: Optional region string (e.g., '22:1000000-2000000') for subsetting
    
    Returns:
        VCF data structure from scikit-allel
    
    Raises:
        FileNotFoundError: If the VCF file doesn't exist
        RuntimeError: If scikit-allel is not installed
    """
    vcf_path = Path(vcf_path)
    
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    
    if not SCIKIT_ALLEL_AVAILABLE:
        raise RuntimeError(
            "scikit-allel not installed. Install via:\n"
            "  pip install scikit-allel"
        )
    
    print(f"Reading VCF file: {vcf_path}")
    
    try:
        if region:
            # Read specific region (if VCF is indexed)
            callset = allel.read_vcf(str(vcf_path), region=region)
        else:
            # Read entire VCF file
            callset = allel.read_vcf(str(vcf_path))
        
        print(f"Successfully loaded VCF file")
        return callset
    except Exception as e:
        raise RuntimeError(f"Error reading VCF file: {e}")


def index_vcf(vcf_path: str) -> None:
    """
    Create a simple index/metadata for a VCF file.
    
    Note: This creates a basic index. For large files, scikit-allel can
    read VCF files directly without needing a tabix index.
    
    Args:
        vcf_path: Path to the VCF file
    """
    vcf_path = Path(vcf_path)
    
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    
    if not SCIKIT_ALLEL_AVAILABLE:
        raise RuntimeError(
            "scikit-allel not installed. Install via:\n"
            "  pip install scikit-allel"
        )
    
    print(f"Reading VCF metadata: {vcf_path}")
    
    try:
        # Read just the header to create metadata
        callset = allel.read_vcf(str(vcf_path), fields=['samples', 'variants/CHROM', 'variants/POS'])
        print(f"VCF file has {len(callset['variants/CHROM'])} variants and {len(callset['samples'])} samples")
        # print("Note: scikit-allel can read VCF files directly without tabix indexing")
    except Exception as e:
        raise RuntimeError(f"Error reading VCF file: {e}")



if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python -m src.utils.vcf_indexer <path_to_vcf> [region]")
        print("  Example: python -m src.utils.vcf_indexer data/file.vcf.gz")
        print("  Example: python -m src.utils.vcf_indexer data/file.vcf.gz 22:1000000-2000000")
        sys.exit(1)
    
    try:
        vcf_path = sys.argv[1]
        region = sys.argv[2] if len(sys.argv) > 2 else None
        
        if region:
            read_vcf(vcf_path, region=region)
        else:
            index_vcf(vcf_path)
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        logger.error("%s", e)
        sys.exit(1)

