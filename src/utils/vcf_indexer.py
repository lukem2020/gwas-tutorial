"""Utility script to work with VCF files using scikit-allel (Python-only, no external tools needed)."""

import sys
from pathlib import Path
import logging
import json
import pickle
from typing import Dict, List, Tuple, Optional

try:
    import yaml

    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False

try:
    import allel
    import numpy as np

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
            "scikit-allel not installed. Install via:\n" "  pip install scikit-allel"
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


def index_vcf(vcf_path: str, index_format: str = "pickle", force: bool = False) -> Path:
    """
    Create a position index for a VCF file to enable fast region queries.

    Creates an index file that maps chromosome positions to variant indices,
    allowing efficient region-based queries without reading the entire VCF.

    Args:
        vcf_path: Path to the VCF file
        index_format: Format for index file ('pickle' or 'json')
        force: Recreate index even if it exists

    Returns:
        Path to the created index file

    Raises:
        FileNotFoundError: If the VCF file doesn't exist
        RuntimeError: If scikit-allel is not installed or indexing fails
    """
    vcf_path = Path(vcf_path)

    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    if not SCIKIT_ALLEL_AVAILABLE:
        raise RuntimeError(
            "scikit-allel not installed. Install via:\n" "  pip install scikit-allel"
        )

    # Determine index file path
    if index_format == "pickle":
        index_path = vcf_path.with_suffix(vcf_path.suffix + ".idx")
    else:
        index_path = vcf_path.with_suffix(vcf_path.suffix + ".idx.json")

    # Check if index already exists
    if index_path.exists() and not force:
        print(f"✓ Index already exists: {index_path}")
        return index_path

    print(f"Creating index for VCF file: {vcf_path}")
    print("This may take a few minutes for large files...")

    try:
        # Read variant positions (this is the slow part for large files)
        print("Reading variant positions...")
        callset = allel.read_vcf(
            str(vcf_path),
            fields=["variants/CHROM", "variants/POS", "variants/ID", "samples"],
        )

        n_variants = len(callset["variants/CHROM"])
        n_samples = len(callset.get("samples", [])) if "samples" in callset else 0

        print(f"Found {n_variants:,} variants")

        # Build position index
        print("Building position index...")
        index_data = {
            "chromosomes": {},
            "variant_count": n_variants,
            "sample_count": n_samples,
            "vcf_path": str(vcf_path),
            "index_format_version": "1.0",
        }

        # Group variants by chromosome
        chroms = callset["variants/CHROM"]
        positions = callset["variants/POS"]
        variant_ids = callset.get("variants/ID", [None] * n_variants)

        # Convert to numpy arrays if needed
        if not isinstance(positions, np.ndarray):
            positions = np.array(positions)
        if not isinstance(chroms, np.ndarray):
            chroms = np.array(chroms)

        # Build index: chromosome -> sorted list of (position, variant_index)
        for i in range(n_variants):
            chrom = str(chroms[i])
            pos = int(positions[i])
            var_id = variant_ids[i] if variant_ids[i] is not None else f"var_{i}"

            if chrom not in index_data["chromosomes"]:
                index_data["chromosomes"][chrom] = []

            index_data["chromosomes"][chrom].append(
                {
                    "position": pos,
                    "variant_index": i,
                    "variant_id": str(var_id) if var_id else None,
                }
            )

        # Sort positions within each chromosome for binary search
        for chrom in index_data["chromosomes"]:
            index_data["chromosomes"][chrom].sort(key=lambda x: x["position"])

        # Save index file
        print(f"Saving index to: {index_path}")
        if index_format == "pickle":
            with open(index_path, "wb") as f:
                pickle.dump(index_data, f, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open(index_path, "w") as f:
                json.dump(index_data, f, indent=2)

        print(f"✓ Index created successfully!")
        print(f"  Variants indexed: {n_variants:,}")
        print(f"  Chromosomes: {len(index_data['chromosomes'])}")
        print(f"  Index file: {index_path}")

        return index_path

    except Exception as e:
        raise RuntimeError(f"Error creating index: {e}")


def load_index(vcf_path: str) -> Optional[Dict]:
    """
    Load an existing index file for a VCF file.

    Args:
        vcf_path: Path to the VCF file

    Returns:
        Index data dictionary or None if index doesn't exist
    """
    vcf_path = Path(vcf_path)

    # Try pickle format first
    index_path = vcf_path.with_suffix(vcf_path.suffix + ".idx")
    if index_path.exists():
        try:
            with open(index_path, "rb") as f:
                return pickle.load(f)
        except Exception:
            pass

    # Try JSON format
    index_path = vcf_path.with_suffix(vcf_path.suffix + ".idx.json")
    if index_path.exists():
        try:
            with open(index_path, "r") as f:
                return json.load(f)
        except Exception:
            pass

    return None


def find_variants_in_region(
    vcf_path: str, region: str, use_index: bool = True
) -> List[int]:
    """
    Find variant indices in a genomic region using the index.

    Args:
        vcf_path: Path to the VCF file
        region: Region string (e.g., '22:1000000-2000000' or '22:1000000')
        use_index: Use index file if available (faster)

    Returns:
        List of variant indices in the region
    """
    vcf_path = Path(vcf_path)

    # Parse region
    if ":" in region:
        parts = region.split(":")
        chrom = parts[0]
        if "-" in parts[1]:
            start, end = map(int, parts[1].split("-"))
        else:
            start = int(parts[1])
            end = start + 1
    else:
        raise ValueError(
            f"Invalid region format: {region}. Use 'chr:start-end' or 'chr:pos'"
        )

    if use_index:
        index_data = load_index(str(vcf_path))
        if index_data:
            chrom_str = str(chrom)
            if chrom_str in index_data["chromosomes"]:
                variants = index_data["chromosomes"][chrom_str]
                # Binary search for variants in range
                indices = [
                    v["variant_index"] for v in variants if start <= v["position"] < end
                ]
                return indices

    # Fallback: read entire VCF (slower)
    callset = allel.read_vcf(str(vcf_path), fields=["variants/CHROM", "variants/POS"])
    chroms = callset["variants/CHROM"]
    positions = callset["variants/POS"]

    indices = [
        i
        for i, (c, p) in enumerate(zip(chroms, positions))
        if str(c) == str(chrom) and start <= int(p) < end
    ]

    return indices


if __name__ == "__main__":
    # Try to load VCF path from config if no arguments provided
    vcf_path = None
    if len(sys.argv) < 2:
        # Try to load from config.yaml
        config_path = Path("config.yaml")
        if config_path.exists() and YAML_AVAILABLE:
            try:
                with open(config_path, "r") as f:
                    config = yaml.safe_load(f)
                if "file_paths" in config and "vcf_file" in config["file_paths"]:
                    vcf_path = config["file_paths"]["vcf_file"]
                    print(f"Using VCF file from config.yaml: {vcf_path}")
            except Exception as e:
                logger.warning(f"Could not load config: {e}")

        if not vcf_path:
            print(
                "Usage: python -m src.utils.vcf_indexer <path_to_vcf> [region] [--force]"
            )
            print("  Example: python -m src.utils.vcf_indexer data/file.vcf.gz")
            print("           (creates index file)")
            print(
                "  Example: python -m src.utils.vcf_indexer data/file.vcf.gz 22:1000000-2000000"
            )
            print("           (reads specific region)")
            print("  Example: python -m src.utils.vcf_indexer data/file.vcf.gz --force")
            print("           (recreates index even if it exists)")
            print("\nOr run without arguments to use VCF file from config.yaml")
            sys.exit(1)
    else:
        vcf_path = sys.argv[1]

    try:
        force = "--force" in sys.argv
        region = None

        # Find region argument (not --force)
        for arg in sys.argv[2:]:
            if arg != "--force" and ":" in arg:
                region = arg
                break

        if region:
            # Read specific region
            print(f"Reading region {region} from {vcf_path}...")
            callset = read_vcf(vcf_path, region=region)
            if callset:
                n_variants = len(callset.get("variants/POS", []))
                print(f"✓ Found {n_variants:,} variants in region {region}")
        else:
            # Create index
            index_path = index_vcf(vcf_path, force=force)
            print(f"\n✓ Indexing complete! Index saved to: {index_path}")
            print("\nYou can now use region queries for faster access:")
            print(f"  python -m src.utils.vcf_indexer {vcf_path} 22:1000000-2000000")

    except (FileNotFoundError, ValueError, RuntimeError) as e:
        logger.error("%s", e)
        sys.exit(1)
