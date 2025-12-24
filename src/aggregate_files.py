"""
Process GEO expression data: Extract probe-to-gene mappings and create gene-level expression matrices.

This module:
1. Parses SOFT files to extract probe-to-gene mappings
2. Aggregates multiple probes per gene (mean expression)
3. Merges expression data with phenotype data
4. Saves analysis-ready gene-level datasets
"""

import pandas as pd
import numpy as np
import gzip
from pathlib import Path
from typing import Dict, Optional, Tuple
import yaml
import warnings
import logging

warnings.filterwarnings("ignore")

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def parse_soft_platform_table(soft_file: str) -> pd.DataFrame:
    """
    Parse the platform table from a SOFT file to extract probe-to-gene mappings.

    Parameters
    ----------
    soft_file : str
        Path to SOFT file (can be .gz compressed)

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: ID (probe), GENE_SYMBOL, GENE (Entrez ID), GENE_NAME, etc.
    """
    soft_file = Path(soft_file)

    # Check if file exists
    if not soft_file.exists():
        logger.error("SOFT file not found: %s", soft_file)
        raise FileNotFoundError("SOFT file not found: %s" % soft_file)

    # Read file (handle both compressed and uncompressed)
    try:
        if soft_file.suffix == ".gz":
            with gzip.open(soft_file, "rt", encoding="utf-8", errors="ignore") as f:
                lines = f.readlines()
        else:
            with open(soft_file, "rt", encoding="utf-8", errors="ignore") as f:
                lines = f.readlines()
    except Exception as e:
        logger.error("Error reading SOFT file %s: %s", soft_file, e)
        raise IOError("Error reading SOFT file: %s" % e)

    # Find platform table section
    table_start = None
    table_end = None
    header_line = None

    for i, line in enumerate(lines):
        if line.startswith("!platform_table_begin"):
            table_start = i + 1
        elif line.startswith("!platform_table_end"):
            table_end = i
            break
        elif table_start is not None and header_line is None:
            # First non-comment line after table_begin is the header
            if line.strip() and not line.startswith("!"):
                header_line = i
                break

    if table_start is None:
        logger.error("Could not find platform table in SOFT file: %s", soft_file)
        raise ValueError("Could not find platform table in SOFT file")

    if header_line is None:
        # Try to find header manually
        for i in range(table_start, min(table_start + 10, len(lines))):
            if lines[i].strip() and not lines[i].startswith("!"):
                header_line = i
                break

    if header_line is None:
        logger.error("Could not find platform table header in SOFT file: %s", soft_file)
        raise ValueError("Could not find platform table header")

    # Parse header
    header = lines[header_line].strip().split("\t")
    header = [col.strip() for col in header]

    # Parse data rows
    data_rows = []
    start_idx = header_line + 1
    end_idx = table_end if table_end else len(lines)

    for i in range(start_idx, end_idx):
        line = lines[i].strip()
        if not line or line.startswith("!"):
            continue
        if line.startswith("!platform_table_end"):
            break

        parts = line.split("\t")
        if len(parts) >= len(header):
            # Pad or truncate to match header length
            row = parts[: len(header)]
            if len(row) < len(header):
                row.extend([""] * (len(header) - len(row)))
            data_rows.append(row)

    if not data_rows:
        logger.error("No data rows found in platform table: %s", soft_file)
        raise ValueError("No data rows found in platform table")

    # Create DataFrame
    platform_df = pd.DataFrame(data_rows, columns=header)

    # Clean up: remove empty rows
    platform_df = platform_df[platform_df.iloc[:, 0].astype(str).str.strip() != ""]

    return platform_df


def extract_probe_to_gene_mapping(platform_df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract probe-to-gene mapping from platform DataFrame.

    Parameters
    ----------
    platform_df : pd.DataFrame
        Platform table DataFrame from parse_soft_platform_table

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: probe_id, gene_symbol, entrez_gene_id, gene_name
    """
    # Determine which columns exist
    probe_col = None
    gene_symbol_col = None
    entrez_col = None
    gene_name_col = None

    # Common column name variations
    for col in platform_df.columns:
        col_upper = col.upper()
        if probe_col is None and ("ID" in col_upper or "PROBE" in col_upper):
            probe_col = col
        if gene_symbol_col is None and (
            "GENE_SYMBOL" in col_upper or "SYMBOL" in col_upper
        ):
            gene_symbol_col = col
        # Check for gene_assignment column (used in some platforms like Illumina)
        if gene_symbol_col is None and "GENE_ASSIGNMENT" in col_upper:
            gene_symbol_col = col
        if entrez_col is None and (
            "GENE" in col_upper
            and "SYMBOL" not in col_upper
            and "NAME" not in col_upper
            and "ASSIGNMENT" not in col_upper
        ):
            entrez_col = col
        if gene_name_col is None and "GENE_NAME" in col_upper:
            gene_name_col = col

    if probe_col is None:
        # Use first column as probe ID
        probe_col = platform_df.columns[0]

    if gene_symbol_col is None:
        logger.error("Could not find GENE_SYMBOL or gene_assignment column in platform table. Available columns: %s", list(platform_df.columns))
        raise ValueError(
            "Could not find GENE_SYMBOL or gene_assignment column in platform table. Available columns: %s" % list(platform_df.columns)
        )

    # Create mapping DataFrame
    # Handle gene_assignment format (e.g., "NM_001005240 // OR4F17 // gene description // location // geneID /// ...")
    if "ASSIGNMENT" in gene_symbol_col.upper():
        # Parse gene_assignment column - format can have multiple entries separated by "///"
        # Each entry: "RefSeq // GeneSymbol // Description // Location // GeneID"
        def extract_gene_symbols(assignment_str):
            """Extract all unique gene symbols from gene_assignment string."""
            if pd.isna(assignment_str) or assignment_str in ["---", "", "N/A", "NA"]:
                return ""
            # Split by "///" to get individual entries
            entries = str(assignment_str).split("///")
            gene_symbols = []
            for entry in entries:
                entry = entry.strip()
                if "//" in entry:
                    parts = [p.strip() for p in entry.split("//")]
                    if len(parts) >= 2 and parts[1]:  # Second part is gene symbol
                        gene_symbols.append(parts[1])
            # Return unique gene symbols joined by "; " if multiple
            unique_symbols = list(
                dict.fromkeys(gene_symbols)
            )  # Preserve order, remove duplicates
            return "; ".join(unique_symbols) if unique_symbols else ""

        gene_symbols_parsed = (
            platform_df[gene_symbol_col].astype(str).apply(extract_gene_symbols)
        )
        mapping_df = pd.DataFrame(
            {
                "probe_id": platform_df[probe_col].astype(str).str.strip(),
                "gene_symbol": gene_symbols_parsed,
            }
        )
    else:
        mapping_df = pd.DataFrame(
            {
                "probe_id": platform_df[probe_col].astype(str).str.strip(),
                "gene_symbol": platform_df[gene_symbol_col].astype(str).str.strip(),
            }
        )

    # Add optional columns
    if entrez_col:
        mapping_df["entrez_gene_id"] = platform_df[entrez_col].astype(str).str.strip()
    else:
        mapping_df["entrez_gene_id"] = ""

    if gene_name_col:
        mapping_df["gene_name"] = platform_df[gene_name_col].astype(str).str.strip()
    else:
        mapping_df["gene_name"] = ""

    # Filter out probes without gene symbols
    # Empty strings, '---', 'N/A', etc.
    invalid_symbols = ["", "---", "N/A", "NA", "nan", "NaN", "NULL", "null"]
    mapping_df = mapping_df[~mapping_df["gene_symbol"].isin(invalid_symbols)]
    mapping_df = mapping_df[mapping_df["gene_symbol"].notna()]
    mapping_df = mapping_df[mapping_df["gene_symbol"].str.strip() != ""]

    # Remove duplicates (keep first)
    mapping_df = mapping_df.drop_duplicates(subset=["probe_id"], keep="first")

    return mapping_df


def aggregate_probes_to_genes(
    expression_df: pd.DataFrame,
    probe_to_gene: pd.DataFrame,
    aggregation_method: str = "mean",
) -> pd.DataFrame:
    """
    Aggregate multiple probes per gene to create gene-level expression matrix.

    Parameters
    ----------
    expression_df : pd.DataFrame
        Expression matrix with probe IDs as index and samples as columns
    probe_to_gene : pd.DataFrame
        Probe-to-gene mapping with columns: probe_id, gene_symbol
    aggregation_method : str
        Method to aggregate multiple probes: 'mean', 'median', or 'max'

    Returns
    -------
    pd.DataFrame
        Gene-level expression matrix with gene symbols as index
    """
    # Reset index if needed
    if expression_df.index.name != "probe_id":
        expression_df = expression_df.reset_index()
        if "ID_REF" in expression_df.columns:
            expression_df = expression_df.rename(columns={"ID_REF": "probe_id"})
        if "probe_id" not in expression_df.columns:
            # Assume first column is probe ID
            first_col = expression_df.columns[0]
            expression_df = expression_df.rename(columns={first_col: "probe_id"})
        expression_df = expression_df.set_index("probe_id")

    # Convert index to string for matching
    expression_df.index = expression_df.index.astype(str).str.strip()
    probe_to_gene["probe_id"] = probe_to_gene["probe_id"].astype(str).str.strip()

    # Handle probes that map to multiple genes (separated by "; ")
    # Expand these so each gene gets the probe's expression value
    probe_to_gene_expanded = []
    for _, row in probe_to_gene[["probe_id", "gene_symbol"]].iterrows():
        probe_id = row["probe_id"]
        gene_symbols = str(row["gene_symbol"]).split("; ")
        for gene_symbol in gene_symbols:
            gene_symbol = gene_symbol.strip()
            if gene_symbol and gene_symbol not in ["", "---", "N/A", "NA"]:
                probe_to_gene_expanded.append(
                    {"probe_id": probe_id, "gene_symbol": gene_symbol}
                )

    probe_to_gene_expanded_df = pd.DataFrame(probe_to_gene_expanded)

    # Merge expression with gene mapping
    expression_with_genes = expression_df.merge(
        probe_to_gene_expanded_df, left_index=True, right_on="probe_id", how="inner"
    )

    # Drop probe_id column and set gene_symbol as index
    expression_with_genes = expression_with_genes.drop(columns=["probe_id"])

    # Group by gene_symbol and aggregate
    if aggregation_method == "mean":
        gene_expression = expression_with_genes.groupby("gene_symbol").mean()
    elif aggregation_method == "median":
        gene_expression = expression_with_genes.groupby("gene_symbol").median()
    elif aggregation_method == "max":
        gene_expression = expression_with_genes.groupby("gene_symbol").max()
    else:
        logger.error("Unknown aggregation method: %s", aggregation_method)
        raise ValueError("Unknown aggregation method: %s" % aggregation_method)

    # Sort by gene symbol
    gene_expression = gene_expression.sort_index()

    return gene_expression


def merge_expression_with_phenotypes(
    gene_expression: pd.DataFrame,
    phenotypes_df: pd.DataFrame,
    sample_id_col: str = "sample_id",
) -> pd.DataFrame:
    """
    Merge gene expression matrix with phenotype data.

    Parameters
    ----------
    gene_expression : pd.DataFrame
        Gene expression matrix (genes × samples)
    phenotypes_df : pd.DataFrame
        Phenotype DataFrame with sample IDs
    sample_id_col : str
        Column name in phenotypes_df that contains sample IDs

    Returns
    -------
    pd.DataFrame
        Combined DataFrame with expression and phenotypes
    """
    # Ensure sample_id column exists
    if sample_id_col not in phenotypes_df.columns:
        logger.error("Column '%s' not found in phenotypes DataFrame", sample_id_col)
        raise ValueError("Column '%s' not found in phenotypes DataFrame" % sample_id_col)

    # Get sample IDs from expression (columns) and phenotypes
    expr_samples = set(gene_expression.columns)
    pheno_samples = set(phenotypes_df[sample_id_col].astype(str))

    # Find common samples
    common_samples = expr_samples.intersection(pheno_samples)

    if len(common_samples) == 0:
        logger.error("No matching samples found between expression and phenotype data")
        raise ValueError(
            "No matching samples found between expression and phenotype data"
        )

    # Filter to common samples
    gene_expression_filtered = gene_expression[list(common_samples)]
    phenotypes_filtered = phenotypes_df[
        phenotypes_df[sample_id_col].astype(str).isin(common_samples)
    ].copy()

    # Transpose expression to have samples as rows
    expression_t = gene_expression_filtered.T
    expression_t.index.name = sample_id_col
    expression_t = expression_t.reset_index()

    # Merge
    combined_df = expression_t.merge(phenotypes_filtered, on=sample_id_col, how="inner")

    return combined_df


def process_geo_dataset(
    geo_accession: str,
    data_dir: str = "data/expression",
    aggregation_method: str = "mean",
) -> Dict:
    """
    Process a single GEO dataset: extract probe-to-gene mappings and create gene-level expression.

    Parameters
    ----------
    geo_accession : str
        GEO accession number (e.g., "GSE28042")
    data_dir : str
        Directory containing GEO data files
    aggregation_method : str
        Method to aggregate probes: 'mean', 'median', or 'max'

    Returns
    -------
    dict
        Dictionary containing:
        - 'probe_to_gene': probe-to-gene mapping DataFrame
        - 'gene_expression': gene-level expression matrix
        - 'combined': merged expression + phenotypes DataFrame
        - 'stats': processing statistics
    """
    data_dir = Path(data_dir)

    # File paths
    soft_file = data_dir / "%s_family.soft.gz" % geo_accession
    expression_file = data_dir / "%s_expression.csv" % geo_accession
    phenotype_file = data_dir / "%s_phenotypes.csv" % geo_accession

    # Check files exist
    missing_files = []
    if not soft_file.exists():
        missing_files.append(soft_file.name)
    if not expression_file.exists():
        missing_files.append(expression_file.name)
    if not phenotype_file.exists():
        missing_files.append(phenotype_file.name)

    if missing_files:
        logger.error("Missing files for %s: %s", geo_accession, ', '.join(missing_files))
        raise FileNotFoundError(
            "Missing files for %s: %s" % (geo_accession, ', '.join(missing_files))
        )

    # Step 1: Parse platform table from SOFT file
    platform_df = parse_soft_platform_table(str(soft_file))

    # Step 2: Extract probe-to-gene mapping
    probe_to_gene = extract_probe_to_gene_mapping(platform_df)

    # Step 3: Load expression data
    try:
        expression_df = pd.read_csv(expression_file, index_col=0)
    except Exception as e:
        logger.error("Error loading expression data from %s: %s", expression_file, e)
        raise

    # Step 4: Aggregate probes to genes
    gene_expression = aggregate_probes_to_genes(
        expression_df, probe_to_gene, aggregation_method
    )

    # Step 5: Load phenotype data
    try:
        phenotypes_df = pd.read_csv(phenotype_file)
    except Exception as e:
        logger.error("Error loading phenotype data from %s: %s", phenotype_file, e)
        raise

    # Step 6: Merge expression with phenotypes
    combined_df = merge_expression_with_phenotypes(gene_expression, phenotypes_df)

    # Calculate statistics
    stats = {
        "n_probes": len(expression_df),
        "n_genes": len(gene_expression),
        "n_samples": len(combined_df),
        "n_phenotype_cols": len(phenotypes_df.columns),
        "probes_per_gene_mean": (
            len(expression_df) / len(gene_expression) if len(gene_expression) > 0 else 0
        ),
    }

    return {
        "probe_to_gene": probe_to_gene,
        "gene_expression": gene_expression,
        "combined": combined_df,
        "stats": stats,
    }


def save_processed_data(
    results: Dict, geo_accession: str, output_dir: str = "data/expression"
) -> Dict[str, Path]:
    """
    Save processed data to files.

    Parameters
    ----------
    results : dict
        Results dictionary from process_geo_dataset
    geo_accession : str
        GEO accession number
    output_dir : str
        Output directory

    Returns
    -------
    dict
        Dictionary mapping file type to file path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    saved_files = {}

    # Save probe-to-gene mapping
    try:
        probe_to_gene_file = output_dir / "%s_probe_to_gene.csv" % geo_accession
        results["probe_to_gene"].to_csv(probe_to_gene_file, index=False)
        saved_files["probe_to_gene"] = probe_to_gene_file
    except Exception as e:
        logger.error("Error saving probe-to-gene mapping for %s: %s", geo_accession, e)
        raise

    # Save gene-level expression matrix
    try:
        gene_expr_file = output_dir / "%s_gene_expression.csv" % geo_accession
        results["gene_expression"].to_csv(gene_expr_file)
        saved_files["gene_expression"] = gene_expr_file
    except Exception as e:
        logger.error("Error saving gene expression for %s: %s", geo_accession, e)
        raise

    # Save combined expression + phenotypes
    try:
        combined_file = output_dir / "%s_gene_expression_with_phenotypes.csv" % geo_accession
        results["combined"].to_csv(combined_file, index=False)
        saved_files["combined"] = combined_file
    except Exception as e:
        logger.error("Error saving combined data for %s: %s", geo_accession, e)
        raise

    return saved_files


def process_all_datasets(
    config_path: str = "config.yaml",
    data_dir: str = "data/expression",
    aggregation_method: str = "mean",
) -> Dict[str, Dict]:
    """
    Process all GEO datasets listed in config.yaml.

    Parameters
    ----------
    config_path : str
        Path to config.yaml
    data_dir : str
        Directory containing GEO data files
    aggregation_method : str
        Method to aggregate probes: 'mean', 'median', or 'max'

    Returns
    -------
    dict
        Dictionary mapping GEO accession to results
    """
    # Load config
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    geo_datasets = config.get("GEO_data_sets", [])

    if not geo_datasets:
        raise ValueError("No GEO datasets found in config.yaml")

    print("=" * 60)
    print("Process GEO Expression Data: Probe-to-Gene Mapping")
    print("=" * 60)
    print("\nFound %d datasets to process:" % len(geo_datasets))
    for gse in geo_datasets:
        print("  - %s" % gse)
    print()

    results = {}

    for i, geo_accession in enumerate(geo_datasets, 1):
        try:
            print("\n[%d/%d] Processing %s..." % (i, len(geo_datasets), geo_accession))
            result = process_geo_dataset(geo_accession, data_dir, aggregation_method)

            # Save files
            saved_files = save_processed_data(result, geo_accession, data_dir)
            result["saved_files"] = saved_files

            results[geo_accession] = result

        except Exception as e:
            logger.error("Error processing %s: %s", geo_accession, e)
            results[geo_accession] = {"error": str(e)}

    return results


def main():
    """Main function to process all GEO datasets from config."""
    try:
        results = process_all_datasets()

        # Print summary
        print("\n" + "=" * 60)
        print("Processing Summary")
        print("=" * 60)

        successful = [gse for gse, data in results.items() if "error" not in data]
        failed = [gse for gse, data in results.items() if "error" in data]

        print("\nSuccessfully processed: %d/%d" % (len(successful), len(results)))
        for gse in successful:
            stats = results[gse]["stats"]
            print("  OK %s: %d genes, %d samples" % (gse, stats['n_genes'], stats['n_samples']))

        if failed:
            print("\nFailed: %d" % len(failed))
            for gse in failed:
                print("  X %s: %s" % (gse, results[gse].get('error', 'Unknown error')))

        print("\n" + "=" * 60)
        print("Output Files")
        print("=" * 60)
        print("For each dataset, the following files are created:")
        print("  - {GSE}_probe_to_gene.csv (probe-to-gene mapping)")
        print("  - {GSE}_gene_expression.csv (gene-level expression matrix)")
        print("  - {GSE}_gene_expression_with_phenotypes.csv (combined dataset)")
        print("\nFiles saved to: data/expression/")

        return 0 if not failed else 1

    except Exception as e:
        logger.error("Error in main: %s", e)
        return 1


if __name__ == "__main__":
    import sys

    sys.exit(main())
