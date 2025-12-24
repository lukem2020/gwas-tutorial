"""
Clean and process GEO expression data: Extract phenotype variables and create binary disease status.

This module:
1. Loads gene expression with phenotypes datasets
2. Parses Characteristics Ch1 to extract phenotype variables
3. Creates binary disease_status column (0 or 1)
4. Removes Characteristics Ch1 column
5. Saves cleaned datasets
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Optional
import yaml
import warnings

warnings.filterwarnings('ignore')


def parse_characteristics(characteristics_str):
    """
    Parse the Characteristics Ch1 string to extract key-value pairs.
    
    Format: "key1: value1; key2: value2; ..."
    
    Parameters
    ----------
    characteristics_str : str
        Characteristics string from GEO data
    
    Returns
    -------
    dict
        Dictionary with parsed key-value pairs
    """
    if pd.isna(characteristics_str) or characteristics_str == '':
        return {}
    
    result = {}
    # Split by semicolon
    pairs = str(characteristics_str).split(';')
    for pair in pairs:
        pair = pair.strip()
        if ':' in pair:
            key, value = pair.split(':', 1)
            key = key.strip().lower()
            value = value.strip()
            result[key] = value
    return result


def extract_phenotype_variables(combined_df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract all phenotype variables from Characteristics Ch1 column.
    
    Parameters
    ----------
    combined_df : pd.DataFrame
        Combined dataset with Characteristics Ch1 column
    
    Returns
    -------
    pd.DataFrame
        DataFrame with extracted phenotype columns added
    """
    print("    Extracting all phenotype variables from Characteristics Ch1...")
    
    # Collect all unique keys across all samples
    all_phenotype_keys = set()
    for idx, row in combined_df.iterrows():
        char_dict = parse_characteristics(row.get('Characteristics Ch1', ''))
        all_phenotype_keys.update(char_dict.keys())
    
    if not all_phenotype_keys:
        print("      No phenotype variables found in Characteristics Ch1")
        return combined_df
    
    print(f"      Found {len(all_phenotype_keys)} unique phenotype variables")
    
    # Extract each phenotype variable into separate columns
    phenotype_data = {}
    for key in all_phenotype_keys:
        phenotype_data[key] = []
    
    for idx, row in combined_df.iterrows():
        char_dict = parse_characteristics(row.get('Characteristics Ch1', ''))
        for key in all_phenotype_keys:
            value = char_dict.get(key, None)
            phenotype_data[key].append(value)
    
    # Add phenotype columns to combined dataframe
    for key, values in phenotype_data.items():
        # Create a clean column name (replace spaces with underscores, make lowercase)
        col_name = key.replace(' ', '_').lower()
        combined_df[col_name] = values
    
    return combined_df, all_phenotype_keys


def create_binary_disease_status(combined_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create binary disease_status column (0 or 1) from phenotype data.
    
    Parameters
    ----------
    combined_df : pd.DataFrame
        Combined dataset with extracted phenotype columns
    
    Returns
    -------
    pd.DataFrame
        DataFrame with binary disease_status column
    """
    print("    Creating binary disease_status column (0 or 1)...")
    
    disease_status_binary = []
    
    for idx, row in combined_df.iterrows():
        # Check if disease_status column exists (from extracted phenotype data)
        disease_status_val = None
        
        # First, check the extracted disease_status column if it exists
        if 'disease_status' in combined_df.columns:
            status_text = str(row['disease_status']).lower() if pd.notna(row['disease_status']) else ''
            if any(term in status_text for term in ['ipf', 'idiopathic pulmonary fibrosis', 'disease', 'case']):
                disease_status_val = 1
            elif any(term in status_text for term in ['control', 'healthy', 'normal']):
                disease_status_val = 0
        
        # If not found, check disease state column
        if disease_status_val is None and 'disease_state' in combined_df.columns:
            status_text = str(row['disease_state']).lower() if pd.notna(row['disease_state']) else ''
            if any(term in status_text for term in ['ipf', 'idiopathic pulmonary fibrosis', 'disease', 'case']):
                disease_status_val = 1
            elif any(term in status_text for term in ['control', 'healthy', 'normal']):
                disease_status_val = 0
        
        # If still not found, check Description or Title columns
        if disease_status_val is None:
            desc = str(row.get('Description', '')).lower()
            title = str(row.get('Title', '')).lower()
            if any(term in desc or term in title for term in ['ipf', 'disease', 'case']):
                disease_status_val = 1
            elif any(term in desc or term in title for term in ['control', 'healthy', 'normal', 'n_']):
                disease_status_val = 0
        
        disease_status_binary.append(disease_status_val)
    
    # Replace or add disease_status column with binary values (0 or 1)
    combined_df['disease_status'] = pd.Series(disease_status_binary, dtype='Int64')  # Int64 allows NaN
    
    return combined_df


def clean_geo_dataset(geo_accession: str, data_dir: str = "data/expression") -> Dict:
    """
    Clean a single GEO dataset: extract phenotypes and create binary disease status.
    
    Parameters
    ----------
    geo_accession : str
        GEO accession number (e.g., "GSE28042")
    data_dir : str
        Directory containing GEO data files
    
    Returns
    -------
    dict
        Dictionary containing cleaned dataframe and statistics
    """
    data_dir = Path(data_dir)
    
    print(f"\nCleaning {geo_accession}...")
    print("=" * 60)
    
    # File path
    combined_file = data_dir / f"{geo_accession}_gene_expression_with_phenotypes.csv"
    
    if not combined_file.exists():
        raise FileNotFoundError(f"File not found: {combined_file}")
    
    # Load combined dataset
    print(f"    Loading: {combined_file.name}")
    combined_df = pd.read_csv(combined_file)
    print(f"      Loaded: {combined_df.shape[0]} samples × {combined_df.shape[1]} columns")
    
    # Check if Characteristics Ch1 exists
    if 'Characteristics Ch1' not in combined_df.columns:
        print("    WARNING: 'Characteristics Ch1' column not found. Dataset may already be cleaned.")
        return {
            'cleaned': combined_df,
            'stats': {
                'n_samples': len(combined_df),
                'n_columns': len(combined_df.columns),
                'phenotype_vars_extracted': 0
            }
        }
    
    # Extract phenotype variables
    combined_df, all_phenotype_keys = extract_phenotype_variables(combined_df)
    
    # Create binary disease_status
    combined_df = create_binary_disease_status(combined_df)
    
    # Remove Characteristics Ch1 column
    if 'Characteristics Ch1' in combined_df.columns:
        combined_df = combined_df.drop(columns=['Characteristics Ch1'])
        print("    Removed 'Characteristics Ch1' column")
    
    # Calculate statistics
    phenotype_cols = [col for col in combined_df.columns if col in [k.replace(' ', '_').lower() for k in all_phenotype_keys]]
    
    stats = {
        'n_samples': len(combined_df),
        'n_columns': len(combined_df.columns),
        'phenotype_vars_extracted': len(all_phenotype_keys),
        'disease_status_1': sum(combined_df['disease_status'] == 1) if 'disease_status' in combined_df.columns else 0,
        'disease_status_0': sum(combined_df['disease_status'] == 0) if 'disease_status' in combined_df.columns else 0,
        'disease_status_missing': sum(combined_df['disease_status'].isna()) if 'disease_status' in combined_df.columns else 0
    }
    
    print(f"\n  OK Cleaning complete!")
    print(f"    - {stats['n_samples']} samples")
    print(f"    - {stats['phenotype_vars_extracted']} phenotype variables extracted")
    print(f"    - Disease (1): {stats['disease_status_1']}, Control (0): {stats['disease_status_0']}")
    
    return {
        'cleaned': combined_df,
        'stats': stats,
        'phenotype_keys': all_phenotype_keys
    }


def save_cleaned_data(results: Dict, geo_accession: str, output_dir: str = "data/expression") -> Path:
    """
    Save cleaned dataset to file.
    
    Parameters
    ----------
    results : dict
        Results dictionary from clean_geo_dataset
    geo_accession : str
        GEO accession number
    output_dir : str
        Output directory
    
    Returns
    -------
    Path
        Path to saved file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save cleaned dataset
    output_file = output_dir / f"{geo_accession}_cleaned.csv"
    results['cleaned'].to_csv(output_file, index=False)
    print(f"    OK Saved: {output_file.name}")
    
    return output_file


def clean_all_datasets(config_path: str = "config.yaml",
                       data_dir: str = "data/expression") -> Dict[str, Dict]:
    """
    Clean all GEO datasets listed in config.yaml.
    
    Parameters
    ----------
    config_path : str
        Path to config.yaml
    data_dir : str
        Directory containing GEO data files
    
    Returns
    -------
    dict
        Dictionary mapping GEO accession to results
    """
    # Load config
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    geo_datasets = config.get('GEO_data_sets', [])
    
    if not geo_datasets:
        raise ValueError("No GEO datasets found in config.yaml")
    
    print("=" * 60)
    print("Clean GEO Expression Data: Extract Phenotypes")
    print("=" * 60)
    print(f"\nFound {len(geo_datasets)} datasets to clean:")
    for gse in geo_datasets:
        print(f"  - {gse}")
    print()
    
    results = {}
    
    for i, geo_accession in enumerate(geo_datasets, 1):
        try:
            print(f"\n[{i}/{len(geo_datasets)}] Processing {geo_accession}...")
            result = clean_geo_dataset(geo_accession, data_dir)
            
            # Save cleaned file
            saved_file = save_cleaned_data(result, geo_accession, data_dir)
            result['saved_file'] = saved_file
            
            results[geo_accession] = result
            
        except Exception as e:
            print(f"  X Error processing {geo_accession}: {e}")
            import traceback
            traceback.print_exc()
            results[geo_accession] = {'error': str(e)}
    
    return results


def main():
    """Main function to clean all GEO datasets from config."""
    try:
        results = clean_all_datasets()
        
        # Print summary
        print("\n" + "=" * 60)
        print("Cleaning Summary")
        print("=" * 60)
        
        successful = [gse for gse, data in results.items() if 'error' not in data]
        failed = [gse for gse, data in results.items() if 'error' in data]
        
        print(f"\nSuccessfully cleaned: {len(successful)}/{len(results)}")
        for gse in successful:
            stats = results[gse]['stats']
            print(f"  OK {gse}: {stats['n_samples']} samples, {stats['phenotype_vars_extracted']} phenotype variables")
            print(f"      Disease (1): {stats['disease_status_1']}, Control (0): {stats['disease_status_0']}")
        
        if failed:
            print(f"\nFailed: {len(failed)}")
            for gse in failed:
                print(f"  X {gse}: {results[gse].get('error', 'Unknown error')}")
        
        print("\n" + "=" * 60)
        print("Output Files")
        print("=" * 60)
        print("For each dataset, the following file is created:")
        print("  - {GSE}_cleaned.csv (cleaned dataset with extracted phenotypes and binary disease_status)")
        print("\nFiles saved to: data/expression/")
        
        return 0 if not failed else 1
        
    except Exception as e:
        print(f"\nX Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())

