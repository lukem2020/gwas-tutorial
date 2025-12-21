"""Utility script to index VCF files using pysam (tabix functionality)."""
import pysam
import sys
import subprocess
from pathlib import Path


def index_vcf(vcf_path: str, use_command: bool = False) -> None:
    """
    Index a VCF file using tabix.
    
    Args:
        vcf_path: Path to the VCF file (can be .vcf or .vcf.gz)
        use_command: If True, use tabix command directly; otherwise use pysam
    """
    vcf_path = Path(vcf_path)
    
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    
    if vcf_path.suffix != '.gz':
        print(f"Note: {vcf_path} is not compressed. Tabix requires bgzip-compressed files.")
        print("Consider compressing with: bgzip <file.vcf> > file.vcf.gz")
        return
    
    print(f"Indexing {vcf_path}...")
    
    if use_command:
        # Try using tabix command directly (available via conda)
        try:
            subprocess.run(['tabix', '-p', 'vcf', str(vcf_path)], check=True)
            print(f"Index created: {vcf_path}.tbi")
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("tabix command not found, falling back to pysam...")
            use_command = False
    
    if not use_command:
        # Use pysam.tabix_index
        pysam.tabix_index(str(vcf_path), preset='vcf', force=True)
        print(f"Index created: {vcf_path}.tbi")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python -m src.utils.vcf_indexer <path_to_vcf.gz>")
        sys.exit(1)
    
    index_vcf(sys.argv[1])

