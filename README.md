# GWAS Tutorial: Gene Expression Analysis Pipeline

A complete Python pipeline for analyzing GEO gene expression data, performing differential expression analysis, and pathway enrichment.

## Overview

This pipeline:
1. Downloads GEO datasets (gene expression + phenotypes)
2. Processes expression data (probe-to-gene mapping, aggregation)
3. Cleans data (extracts phenotypes, creates binary disease status)
4. Performs differential expression analysis
5. Visualizes results (volcano plots, heatmaps)
6. Prepares data for pathway enrichment

## Prerequisites

- Python 3.8 or higher
- pip package manager
- Git (optional, for cloning)

## Quick Start

### Step 1: Create Virtual Environment

```bash
# Create virtual environment
python -m venv .venv

# Activate virtual environment
# On Windows:
.venv\Scripts\activate

# On macOS/Linux:
source .venv/bin/activate
```

### Step 2: Install Requirements

```bash
pip install -r requirements.txt
```

### Step 3: Run Data Processing Pipeline

The main pipeline downloads, processes, and cleans the GEO datasets:

```bash
python main.py
```

This will:
- Download GEO datasets listed in `config.yaml`
- Extract probe-to-gene mappings
- Aggregate expression to gene level
- Clean phenotype data
- Save processed files to `data/expression/`

**Expected output files:**
- `{GSE}_expression.csv` - Raw expression data
- `{GSE}_phenotypes.csv` - Phenotype data
- `{GSE}_probe_to_gene.csv` - Probe-to-gene mapping
- `{GSE}_gene_expression.csv` - Gene-level expression matrix
- `{GSE}_gene_expression_with_phenotypes.csv` - Combined dataset
- `{GSE}_cleaned.csv` - Cleaned dataset with extracted phenotypes

### Step 4: Run Analysis Notebook

Open and run the analysis notebook:

```bash
# If using Jupyter Lab:
jupyter lab notebooks/analysis.ipynb

# If using Jupyter Notebook:
jupyter notebook notebooks/analysis.ipynb

# If using VS Code:
# Just open notebooks/analysis.ipynb in VS Code
```

**In the notebook, run cells sequentially:**

1. **Cell 1**: Load data (loads cleaned dataset)
2. **Cell 2**: Pipeline overview (markdown)
3. **Cell 3**: Data preparation and QC
4. **Cell 4**: Differential expression analysis
5. **Cell 5**: Volcano plot visualization
6. **Cell 6**: Heatmap visualization
7. **Cell 7**: Pathway enrichment preparation
8. **Cell 8**: Results export

**Results will be saved to:**
- `results/differential_expression_results.csv` - Full results
- `results/significant_genes_results.csv` - Significant genes only
- `results/significant_genes.txt` - Gene list for enrichment tools
- `results/ranked_genes_for_gsea.txt` - Ranked list for GSEA
- `results/analysis_summary.csv` - Summary statistics

## Project Structure

```
gwas-tutorial/
├── main.py                 # Main pipeline entry point
├── config.yaml             # Configuration (GEO datasets to process)
├── requirements.txt        # Python dependencies
├── README.md              # This file
├── src/                   # Source code modules
│   ├── download_geo_data.py    # Download GEO datasets
│   ├── aggregate_files.py      # Process expression data
│   └── clean_data.py           # Clean phenotype data
├── notebooks/             # Jupyter notebooks
│   └── analysis.ipynb     # Differential expression analysis
├── data/                  # Data directory (created automatically)
│   └── expression/        # Processed expression data
└── results/               # Analysis results (created by notebook)
```

## Configuration

Edit `config.yaml` to specify which GEO datasets to process:

```yaml
GEO_data_sets:
  - GSE28042
  - GSE38958
  - GSE33566
  - GSE93606
```

## Pipeline Steps

### 1. Download (`src/download_geo_data.py`)
- Downloads Series Matrix, SOFT, and MINiML files from GEO
- Extracts expression and phenotype data
- Saves to `data/expression/`

### 2. Process (`src/aggregate_files.py`)
- Parses SOFT files for probe-to-gene mappings
- Aggregates probe-level to gene-level expression
- Merges expression with phenotypes
- Saves gene-level datasets

### 3. Clean (`src/clean_data.py`)
- Extracts phenotype variables from structured text
- Creates binary disease_status column (0/1)
- Removes original Characteristics Ch1 column
- Saves cleaned datasets

### 4. Analyze (`notebooks/analysis.ipynb`)
- Differential expression analysis (disease vs control)
- Statistical testing with FDR correction
- Visualization (volcano plots, heatmaps)
- Pathway enrichment preparation

## Output Files

### From `main.py`:
- `data/expression/{GSE}_expression.csv` - Raw expression
- `data/expression/{GSE}_phenotypes.csv` - Phenotypes
- `data/expression/{GSE}_probe_to_gene.csv` - Probe mapping
- `data/expression/{GSE}_gene_expression.csv` - Gene expression
- `data/expression/{GSE}_gene_expression_with_phenotypes.csv` - Combined
- `data/expression/{GSE}_cleaned.csv` - Cleaned data

### From `analysis.ipynb`:
- `results/differential_expression_results.csv` - All genes
- `results/significant_genes_results.csv` - Significant only
- `results/significant_genes.txt` - For Enrichr/DAVID
- `results/ranked_genes_for_gsea.txt` - For GSEA
- `results/analysis_summary.csv` - Summary stats

## Troubleshooting

### Virtual Environment Issues

**Windows:**
```bash
# If activation doesn't work:
.venv\Scripts\python.exe -m pip install -r requirements.txt
```

**macOS/Linux:**
```bash
# If activation doesn't work:
source .venv/bin/activate
```

### Import Errors

If you get import errors, make sure:
1. Virtual environment is activated
2. All requirements are installed: `pip install -r requirements.txt`
3. You're running from the project root directory

### Data Loading Issues

If notebooks can't find data:
1. Make sure `main.py` has completed successfully
2. Check that files exist in `data/expression/`
3. Verify the dataset name matches what's in the notebook

### Memory Issues

For large datasets:
- Process one dataset at a time
- Close other applications
- Consider using a machine with more RAM

## Next Steps After Analysis

1. **Pathway Enrichment**: Use saved gene lists with:
   - [Enrichr](https://maayanlab.cloud/Enrichr/)
   - [DAVID](https://david.ncifcrf.gov/)
   - [GSEA](https://www.gsea-msigdb.org/gsea/index.jsp)

2. **Literature Review**: Investigate top differentially expressed genes

3. **Integration**: Combine with published GWAS/eQTL results

4. **Validation**: Test findings in independent datasets

## Dependencies

See `requirements.txt` for full list. Key packages:
- `pandas` - Data manipulation
- `numpy` - Numerical computing
- `scipy` - Statistical functions
- `statsmodels` - Multiple testing correction
- `matplotlib` / `seaborn` - Visualization
- `scikit-learn` - Machine learning utilities
- `PyYAML` - Configuration file parsing

4. Verify data files are present



