# GWAS Tutorial: Step-by-Step Pipeline

A complete Python pipeline for GWAS, polygenic risk scoring (PRS), and eQTL analysis using **1000 Genomes Project** genotype data and **GEO** expression data.

## Quick Start

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Installation

1. Clone or download this repository
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

### Running the Pipeline

#### Step 1: Configure Project

Set up your project configuration:

```bash
python step1_configure.py
```

This will:
- Load and validate `config.yaml`
- Create necessary directory structure
- Display configuration summary

The `config.yaml` file contains all settings for:
- Data directories
- 1000 Genomes download settings
- GEO expression data settings
- Quality control thresholds
- GWAS parameters
- PRS and eQTL settings

#### Step 2: Download Genotype Data

Download 1000 Genomes VCF files for specified chromosomes:

```bash
python step2_download_genotypes.py
```

This will:
- Download VCF files from 1000 Genomes FTP server
- Show progress for each chromosome
- Verify downloaded files
- Skip already downloaded files (unless forced)

**Note**: Downloads can be large (100s of MB to several GB per chromosome). 
The script downloads chromosomes specified in `config.yaml` (default: 2, 3, 4, 6, 8, 9, 10, 11).

## Project Structure

```
gwas-tutorial/
├── config.yaml              # Project configuration
├── step1_configure.py       # Step 1: Configure project
├── src/                     # Python modules
│   ├── __init__.py
│   └── config.py           # Configuration management
├── notebooks/              # Jupyter notebooks
│   └── visualization.ipynb
├── data/                   # Data directory (created automatically)
│   ├── genotypes/         # 1000 Genomes VCF files
│   ├── phenotypes/        # Phenotype CSV files
│   └── expression/         # GEO expression data
├── results/                # Analysis results (created automatically)
│   └── plots/             # Generated plots
└── requirements.txt        # Python dependencies
```

## Pipeline Steps

### Pre-Analysis: Data Preparation & Quality Control

1. ✅ **Configure project** - Set up configuration and directories
2. Download genotype data (1000 Genomes)
3. Prepare phenotype data
4. Download expression data (GEO)
5. Load genotypes
6. Load phenotype
7. Quality Control (QC)
8. Calculate PCA

### Post-Analysis: Association Testing & Interpretation

9. Run GWAS
10. Calculate PRS
11. Run eQTL analysis
12. Visualize results
13. Save results
14. Review & interpret

## Configuration

Edit `config.yaml` to customize:
- Chromosomes to analyze
- QC thresholds
- GWAS model parameters
- Output settings

See `guide.md` for detailed best practices and guidelines.

## Data Sources

- **Genotype Data**: [1000 Genomes Project](https://www.internationalgenome.org/)
- **Expression Data**: [GEO (Gene Expression Omnibus)](https://www.ncbi.nlm.nih.gov/geo/)

## License

This tutorial is for educational purposes.

