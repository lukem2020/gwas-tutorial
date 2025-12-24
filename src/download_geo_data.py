"""
Download GEO datasets listed in config.yaml.
Downloads Series Matrix files and extracts expression and phenotype data.
"""

import pandas as pd
import urllib.request
import gzip
from pathlib import Path
from typing import List, Dict, Optional
import yaml
import time
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def download_geo_file(
    geo_accession: str, file_type: str, output_dir: str
) -> Optional[str]:
    """
    Download a GEO file (Series Matrix, SOFT, or MINiML).

    Parameters
    ----------
    geo_accession : str
        GEO accession number (e.g., "GSE12345")
    file_type : str
        Type of file: 'matrix', 'soft', or 'miniml'
    output_dir : str
        Directory to save downloaded file

    Returns
    -------
    str or None
        Path to downloaded file, or None if download failed
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Construct file names and URLs
    first_two = geo_accession[3:5]  # e.g., "93" from "GSE93606"
    base_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE%snnn/%s" % (first_two, geo_accession)

    if file_type == "matrix":
        filename = "%s_series_matrix.txt.gz" % geo_accession
        url = "%s/matrix/%s" % (base_url, filename)
        output_path = output_dir / filename
    elif file_type == "soft":
        filename = "%s_family.soft.gz" % geo_accession
        url = "%s/soft/%s" % (base_url, filename)
        output_path = output_dir / filename
    elif file_type == "miniml":
        filename = "%s_family.xml.tgz" % geo_accession
        url = "%s/miniml/%s" % (base_url, filename)
        output_path = output_dir / filename
    else:
        logger.error("Unknown file_type: %s. Must be 'matrix', 'soft', or 'miniml'", file_type)
        raise ValueError(
            "Unknown file_type: %s. Must be 'matrix', 'soft', or 'miniml'" % file_type
        )

    # Skip if already downloaded
    if output_path.exists():
        file_size = output_path.stat().st_size
        if file_size > 0:
            return str(output_path)

    try:
        urllib.request.urlretrieve(url, output_path)
        file_size = output_path.stat().st_size
        return str(output_path)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            # Try alternative URL format (direct download)
            url_alt = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=%s&format=%s&file=%s" % (geo_accession, file_type, filename)
            try:
                urllib.request.urlretrieve(url_alt, output_path)
                file_size = output_path.stat().st_size
                if file_size > 0:
                    return str(output_path)
            except Exception as e2:
                logger.error("Alternative URL download failed for %s: %s", filename, e2)

            logger.error("File not available (404): %s", filename)
            return None
        else:
            logger.error("HTTP Error %s for %s: %s", e.code, filename, e)
            return None
    except Exception as e:
        logger.error("Download failed for %s: %s", filename, e)
        return None


def download_geo_series_matrix(geo_accession: str, output_dir: str) -> Optional[str]:
    """
    Download GEO Series Matrix file (backward compatibility).

    Parameters
    ----------
    geo_accession : str
        GEO accession number (e.g., "GSE12345")
    output_dir : str
        Directory to save downloaded file

    Returns
    -------
    str or None
        Path to downloaded file, or None if download failed
    """
    return download_geo_file(geo_accession, "matrix", output_dir)


def parse_geo_series_matrix(matrix_file: str) -> Dict:
    """
    Parse GEO Series Matrix file to extract expression and metadata.

    Parameters
    ----------
    matrix_file : str
        Path to Series Matrix file

    Returns
    -------
    dict
        Dictionary with 'expression' (DataFrame) and 'metadata' (DataFrame)
    """
    matrix_file = Path(matrix_file)

    # Read the file
    try:
        with gzip.open(matrix_file, "rt", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
    except Exception as e:
        # Try as plain text
        try:
            with open(matrix_file, "rt", encoding="utf-8", errors="ignore") as f:
                lines = f.readlines()
        except Exception as e2:
            logger.error("Error reading Series Matrix file %s: %s", matrix_file, e2)
            raise

    # Find data section
    data_start = None
    data_end = None
    metadata_lines = []

    for i, line in enumerate(lines):
        if line.startswith("!Series_matrix_table_begin"):
            data_start = i + 1
        elif line.startswith("!Series_matrix_table_end"):
            data_end = i
            break
        elif line.startswith("!"):
            metadata_lines.append(line.strip())

    # Parse metadata
    metadata_dict = {}
    for line in metadata_lines:
        if "=" in line:
            key, value = line[1:].split("=", 1)
            key = key.strip()
            value = value.strip().strip('"')
            if key not in metadata_dict:
                metadata_dict[key] = []
            metadata_dict[key].append(value)

    # Parse expression data
    if data_start is None:
        # Try to find header line with ID_REF
        for i, line in enumerate(lines):
            if "ID_REF" in line or (line.strip() and not line.startswith("!")):
                data_start = i
                break

    if data_start is None:
        logger.error("Could not find data section in Series Matrix file: %s", matrix_file)
        raise ValueError("Could not find data section")

    # Find header and data
    header_line = None
    for i in range(data_start, len(lines)):
        line = lines[i].strip()
        if line and not line.startswith("!"):
            header_line = i
            break

    if header_line is None:
        logger.error("Could not find header line in Series Matrix file: %s", matrix_file)
        raise ValueError("Could not find header line")

    # Read data lines
    end_line = data_end if data_end else len(lines)
    data_lines = []
    for i in range(header_line, end_line):
        line = lines[i].strip()
        if line and not line.startswith("!"):
            data_lines.append(line)

    if not data_lines:
        logger.error("No data lines found in Series Matrix file: %s", matrix_file)
        raise ValueError("No data lines found")

    # Parse header
    header_raw = data_lines[0].split("\t")
    header = [col.strip().strip('"') for col in header_raw]

    # Parse data rows
    data_rows = []
    for line in data_lines[1:]:
        parts = line.split("\t")
        parts = [p.strip().strip('"') for p in parts]
        if len(parts) >= len(header):
            data_rows.append(parts[: len(header)])
        elif len(parts) > 0:
            # Pad if needed
            parts.extend([""] * (len(header) - len(parts)))
            data_rows.append(parts)

    # Create expression DataFrame
    if data_rows:
        expression_df = pd.DataFrame(data_rows, columns=header)
        if len(expression_df.columns) > 0:
            first_col = expression_df.columns[0]
            expression_df = expression_df.set_index(first_col)
        # Convert to numeric
        for col in expression_df.columns:
            expression_df[col] = pd.to_numeric(expression_df[col], errors="coerce")
    else:
        # Empty data - just extract sample IDs from header
        sample_ids = header[1:] if len(header) > 1 else []
        expression_df = pd.DataFrame(index=[], columns=sample_ids)

    # Create metadata DataFrame
    sample_ids = (
        expression_df.columns.tolist() if len(expression_df.columns) > 0 else []
    )
    if not sample_ids and "Series_sample_id" in metadata_dict:
        sample_ids = metadata_dict["Series_sample_id"]

    metadata_df = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "geo_accession": (
                [metadata_dict.get("Series_geo_accession", [""])[0]] * len(sample_ids)
                if sample_ids
                else []
            ),
        }
    )

    return {
        "expression": expression_df,
        "metadata": metadata_df,
        "raw_metadata": metadata_dict,
    }


def extract_phenotypes(geo_data: Dict) -> pd.DataFrame:
    """
    Extract phenotype information from GEO metadata.

    Parameters
    ----------
    geo_data : dict
        Dictionary from parse_geo_series_matrix

    Returns
    -------
    pd.DataFrame
        Phenotype DataFrame
    """
    metadata = geo_data["metadata"]
    raw_metadata = geo_data["raw_metadata"]

    phenotype_df = pd.DataFrame({"sample_id": metadata["sample_id"]})

    # Look for sample characteristics
    found_cols = []
    for key, values in raw_metadata.items():
        key_lower = key.lower()
        if "sample" in key_lower and (
            "characteristic" in key_lower or "title" in key_lower
        ):
            for value in values:
                if ":" in value:
                    parts = value.split(":", 1)
                    if len(parts) == 2:
                        sample_id = parts[0].strip()
                        pheno_value = parts[1].strip()
                        col_name = (
                            key.replace("Series_sample_characteristic_", "")
                            .replace("Series_sample_", "")
                            .replace("_", " ")
                            .title()
                        )
                        if col_name and col_name not in phenotype_df.columns:
                            phenotype_df[col_name] = None
                            found_cols.append(col_name)
                        if col_name:
                            idx = phenotype_df[
                                phenotype_df["sample_id"] == sample_id
                            ].index
                            if len(idx) > 0:
                                phenotype_df.loc[idx[0], col_name] = pheno_value

    if not found_cols:
        phenotype_df["phenotype"] = None

    return phenotype_df


def parse_soft_file(soft_file: str) -> Dict:
    """
    Parse GEO SOFT file to extract detailed metadata and phenotypes.
    SOFT files contain more detailed sample information than Series Matrix.

    Parameters
    ----------
    soft_file : str
        Path to SOFT file

    Returns
    -------
    dict
        Dictionary with parsed metadata
    """
    import gzip

    soft_file = Path(soft_file)

    # Read SOFT file
    try:
        with gzip.open(soft_file, "rt", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
    except Exception as e:
        try:
            with open(soft_file, "rt", encoding="utf-8", errors="ignore") as f:
                lines = f.readlines()
        except Exception as e2:
            logger.error("Error reading SOFT file %s: %s", soft_file, e2)
            raise

    # Parse SOFT format
    # SOFT files have sections: ^SERIES, ^SAMPLE, ^DATABASE
    samples = {}
    current_sample = None
    series_metadata = {}

    for line in lines:
        line = line.strip()
        if line.startswith("^SAMPLE"):
            # New sample section
            parts = line.split("=")
            if len(parts) > 1:
                current_sample = parts[1].strip()
                samples[current_sample] = {}
        elif line.startswith("!Sample_") and current_sample:
            # Sample metadata
            if "=" in line:
                key, value = line[1:].split("=", 1)
                key = key.strip()
                value = value.strip().strip('"')
                if key not in samples[current_sample]:
                    samples[current_sample][key] = []
                samples[current_sample][key].append(value)
        elif line.startswith("!Series_"):
            # Series metadata
            if "=" in line:
                key, value = line[1:].split("=", 1)
                key = key.strip()
                value = value.strip().strip('"')
                if key not in series_metadata:
                    series_metadata[key] = []
                series_metadata[key].append(value)

    return {"samples": samples, "series_metadata": series_metadata}


def download_geo_dataset(geo_accession: str, output_dir: str) -> Dict:
    """
    Download and parse a single GEO dataset.
    Downloads Series Matrix, SOFT, and MINiML files.

    Parameters
    ----------
    geo_accession : str
        GEO accession number
    output_dir : str
        Output directory

    Returns
    -------
    dict
        Dictionary with 'expression', 'metadata', 'phenotypes' DataFrames
    """
    output_dir = Path(output_dir)

    # Download all file types
    matrix_file = download_geo_file(geo_accession, "matrix", str(output_dir))
    soft_file = download_geo_file(geo_accession, "soft", str(output_dir))
    miniml_file = download_geo_file(geo_accession, "miniml", str(output_dir))

    if not matrix_file:
        logger.error("Failed to download Series Matrix for %s", geo_accession)
        raise ValueError("Failed to download Series Matrix for %s" % geo_accession)

    # Parse Series Matrix (required for expression data)
    geo_data = parse_geo_series_matrix(matrix_file)

    # Parse SOFT file if available (has more detailed metadata)
    soft_data = None
    if soft_file:
        try:
            soft_data = parse_soft_file(soft_file)
            geo_data["soft_data"] = soft_data

            # Try to extract better phenotype information from SOFT
            if soft_data and "samples" in soft_data:
                # Update phenotypes with SOFT metadata
                phenotypes = extract_phenotypes_from_soft(geo_data, soft_data)
                geo_data["phenotypes"] = phenotypes
            else:
                # Fall back to Series Matrix extraction
                phenotypes = extract_phenotypes(geo_data)
                geo_data["phenotypes"] = phenotypes
        except Exception as e:
            logger.warning("Could not parse SOFT file for %s: %s", geo_accession, e)
            # Fall back to Series Matrix extraction
            phenotypes = extract_phenotypes(geo_data)
            geo_data["phenotypes"] = phenotypes
    else:
        # No SOFT file, use Series Matrix extraction
        phenotypes = extract_phenotypes(geo_data)
        geo_data["phenotypes"] = phenotypes

    # Save files
    try:
        expression_file = output_dir / "%s_expression.csv" % geo_accession
        phenotype_file = output_dir / "%s_phenotypes.csv" % geo_accession

        geo_data["expression"].to_csv(expression_file)
        phenotypes.to_csv(phenotype_file, index=False)
    except Exception as e:
        logger.error("Error saving files for %s: %s", geo_accession, e)
        raise

    return geo_data


def extract_phenotypes_from_soft(geo_data: Dict, soft_data: Dict) -> pd.DataFrame:
    """
    Extract phenotype information from SOFT file data.
    SOFT files have more detailed sample characteristics.

    Parameters
    ----------
    geo_data : dict
        Dictionary from parse_geo_series_matrix
    soft_data : dict
        Dictionary from parse_soft_file

    Returns
    -------
    pd.DataFrame
        Phenotype DataFrame with extracted information
    """
    metadata = geo_data["metadata"]
    samples = soft_data.get("samples", {})

    # Create phenotype DataFrame with sample IDs
    phenotype_df = pd.DataFrame({"sample_id": metadata["sample_id"]})

    # Extract all sample characteristics from SOFT file
    found_cols = []

    # Get all unique keys across all samples
    all_keys = set()
    for sample_id, sample_data in samples.items():
        all_keys.update(sample_data.keys())

    # Create columns for each characteristic
    for key in all_keys:
        if key.startswith("Sample_"):
            # Clean column name
            col_name = key.replace("Sample_", "").replace("_", " ").title()
            if col_name and col_name not in phenotype_df.columns:
                phenotype_df[col_name] = None
                found_cols.append(col_name)

    # Fill in values
    for idx, sample_id in enumerate(phenotype_df["sample_id"]):
        if sample_id in samples:
            sample_data = samples[sample_id]
            for key, values in sample_data.items():
                if key.startswith("Sample_"):
                    col_name = key.replace("Sample_", "").replace("_", " ").title()
                    if col_name in phenotype_df.columns:
                        # Join multiple values if present
                        value = (
                            "; ".join(str(v) for v in values)
                            if isinstance(values, list)
                            else str(values)
                        )
                        phenotype_df.loc[idx, col_name] = value

    if not found_cols:
        phenotype_df["phenotype"] = None

    return phenotype_df


def download_all_geo_datasets(
    config_path: str = "config.yaml", output_dir: str = "data/expression"
) -> Dict[str, Dict]:
    """
    Download all GEO datasets listed in config.yaml.

    Parameters
    ----------
    config_path : str
        Path to config.yaml file
    output_dir : str
        Directory to save downloaded files

    Returns
    -------
    dict
        Dictionary mapping GEO accession to data dictionaries
    """
    # Load config
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    geo_datasets = config.get("GEO_data_sets", [])

    if not geo_datasets:
        logger.error("No GEO datasets found in config.yaml")
        raise ValueError("No GEO datasets found in config.yaml")

    results = {}
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    for i, geo_accession in enumerate(geo_datasets, 1):
        try:
            geo_data = download_geo_dataset(geo_accession, str(output_path))
            results[geo_accession] = geo_data

            # Small delay to avoid overwhelming the server
            time.sleep(1)

        except Exception as e:
            logger.error("Error processing %s: %s", geo_accession, e)
            results[geo_accession] = {"error": str(e)}

    return results


def main():
    """Main function to download all GEO datasets from config."""
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
        successful = [gse for gse, data in results.items() if "error" not in data]
        failed = [gse for gse, data in results.items() if "error" in data]

        print("Successfully downloaded: %d/%d" % (len(successful), len(results)))
        if successful:
            for gse in successful:
                data = results[gse]
                expr_shape = data["expression"].shape
                n_samples = len(data["phenotypes"])
                print("  OK %s: %d genes × %d samples, %d phenotypes" % (gse, expr_shape[0], expr_shape[1], n_samples))

        if failed:
            print("\nFailed: %d" % len(failed))
            for gse in failed:
                print("  X %s: %s" % (gse, results[gse].get('error', 'Unknown error')))

        print()
        print("Files saved to: data/expression/")
        print("  - {GSE}_expression.csv (expression data)")
        print("  - {GSE}_phenotypes.csv (phenotype data)")

        return 0 if not failed else 1

    except Exception as e:
        logger.error("Error in main: %s", e)
        return 1


if __name__ == "__main__":
    import sys

    sys.exit(main())
