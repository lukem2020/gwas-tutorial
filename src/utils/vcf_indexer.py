"""Utility script to index VCF files using tabix command-line tool."""
import sys
import subprocess
from pathlib import Path
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())


def index_vcf(vcf_path: str) -> None:
    """
    Index a VCF file using tabix command-line tool.
    
    Args:
        vcf_path: Path to the VCF file (must be .vcf.gz compressed format)
    
    Raises:
        FileNotFoundError: If the VCF file doesn't exist
        RuntimeError: If tabix command is not found or indexing fails
    """
    vcf_path = Path(vcf_path)
    
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    
    if vcf_path.suffix != '.gz':
        raise ValueError(
            f"{vcf_path} is not compressed. Tabix requires bgzip-compressed files.\n"
            "Consider compressing with: bgzip <file.vcf> > file.vcf.gz"
        )
    
    print(f"Indexing {vcf_path}...")
    
    try:
        subprocess.run(['tabix', '-p', 'vcf', str(vcf_path)], check=True)
        print(f"Index created: {vcf_path}.tbi")
    except FileNotFoundError:
        raise RuntimeError(
            "tabix command not found. Install via:\n"
            "  conda install -c conda-forge -c bioconda tabix"
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error indexing VCF with tabix: {e}")



if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python -m src.utils.vcf_indexer <path_to_vcf.gz>")
        sys.exit(1)
    
    try:
        index_vcf(sys.argv[1])
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        logger.error("%s", e)
        sys.exit(1)

