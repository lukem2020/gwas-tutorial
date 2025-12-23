#v 1.0
# GWAS + PRS + eQTL Pipeline Guide
## Using 1000 Genomes & GEO Data

This guide provides step-by-step instructions for performing a complete GWAS, Polygenic Risk Scoring (PRS), and eQTL analysis pipeline using genotype data from the 1000 Genomes Project and expression/phenotype data from GEO (Gene Expression Omnibus).

---

## Table of Contents

1. [GEO Data Acquisition](#geo-data-acquisition)
2. [Data Exploration](#data-exploration)
3. [Genotype Data Preparation](#genotype-data-preparation)
4. [Quality Control](#quality-control)
5. [Association Testing](#association-testing)
6. [Polygenic Risk Scoring](#polygenic-risk-scoring)
7. [eQTL Analysis](#eqtl-analysis)

---

## GEO Data Acquisition

### Overview

GEO (Gene Expression Omnibus) is a public repository of gene expression and phenotype data. For this pipeline, we download GEO datasets that contain both expression data and phenotype information (e.g., disease status, clinical variables).

### Step 1: Configure GEO Datasets

Edit `config.yaml` to specify which GEO datasets you want to download:

```yaml
GEO_data_sets:
  - GSE93606 # Host-Microbial interactions in Idiopathic Pulmonary Fibrosis
  - GSE38958 # Profiling of Gene Expression in Idiopathic Pulmonary Fibrosis
  - GSE28042 # Peripheral Blood Mononuclear Cell Gene Expression Profiles
  - GSE33566 # The Peripheral Blood Transcriptome Predicts Disease in IPF
```

**Finding GEO Datasets:**
- Search GEO: https://www.ncbi.nlm.nih.gov/geo/
- Look for datasets with:
  - Expression data (microarray or RNA-seq)
  - Phenotype information (disease status, clinical variables)
  - Sample sizes appropriate for your analysis
  - Relevant tissue types for your research question

### Step 2: Download GEO Datasets

Run the download script to fetch all datasets specified in `config.yaml`:

```bash
python download_geo_datasets.py
```

**What Gets Downloaded:**

For each GEO dataset (e.g., `GSE93606`), the script downloads three file types:

1. **Series Matrix File** (`{GSE}_series_matrix.txt.gz`)
   - Pre-processed expression data matrix
   - Format: Genes/probes × Samples
   - Used for quick access to expression values

2. **SOFT Formatted Family File** (`{GSE}_family.soft.gz`)
   - Complete metadata in SOFT format
   - Contains platform annotation (probe-to-gene mappings)
   - Detailed sample characteristics and phenotypes
   - More comprehensive than Series Matrix

3. **MINiML Formatted Family File** (`{GSE}_family.xml.tgz`)
   - XML format metadata
   - Alternative format for programmatic access
   - Contains same information as SOFT file

**Output Files Created:**

After downloading and parsing, the following files are created in `data/expression/`:

- `{GSE}_expression.csv` - Expression data matrix (genes/probes × samples)
- `{GSE}_phenotypes.csv` - Extracted phenotype metadata
- `{GSE}_series_matrix.txt.gz` - Original Series Matrix file
- `{GSE}_family.soft.gz` - Original SOFT file
- `{GSE}_family.xml.tgz` - Original MINiML file

### Step 3: Understanding the Data Structure

#### Expression Data (`{GSE}_expression.csv`)

- **Format**: CSV file with probe/gene IDs as index and sample IDs (GSM numbers) as columns
- **Values**: Log2-transformed or normalized expression values
- **Example**:
  ```
  ID_REF,GSM2458563,GSM2458564,GSM2458565,...
  7892501,1.268,1.604,1.853,...
  7892502,3.255,1.910,2.784,...
  ```

#### Phenotype Data (`{GSE}_phenotypes.csv`)

- **Format**: CSV file with one row per sample
- **Columns**: 
  - `sample_id` - GEO sample ID (GSM number)
  - `Characteristics Ch1` - Structured phenotype information (semicolon-separated)
  - `Description` - Sample description
  - `Title` - Sample title
  - Additional metadata columns (protocols, contact info, etc.)

- **Phenotype Information**:
  The `Characteristics Ch1` column contains structured data like:
  ```
  tissue: whole blood; disease state: control; gender: male; 
  survival (months): 34; age (years): 71; fvc: NA; dlco: NA; 
  composite_end_point: 0
  ```

### Step 4: Mapping Probe IDs to Gene Symbols

Microarray platforms use probe IDs (e.g., "7892501") that need to be mapped to gene symbols (e.g., "ACTB"). The SOFT file contains platform annotation information.

**Process:**
1. Parse the SOFT file's platform table section
2. Extract probe-to-gene mappings from the `gene_assignment` column
3. Map format: `NM_001005240 // OR4F17 // gene description // ...`
   - Gene symbol is the second field after splitting by `//`
4. Create expression matrix with gene symbols as index
5. Aggregate multiple probes per gene (take mean if multiple probes map to same gene)

**Result:**
- Expression matrix with gene symbols instead of probe IDs
- Enables gene-level analysis and interpretation

### Step 5: Data Exploration

Use the Jupyter notebook `notebooks/data_exploration.ipynb` to:

1. **Explore Phenotype Data:**
   - View all columns and data types
   - Parse structured phenotype information
   - Identify controls vs. cases
   - Check for missing values

2. **Explore Expression Data:**
   - View expression matrix structure
   - Identify highly expressed genes
   - Check expression value distributions
   - Verify sample matching between expression and phenotype data

3. **Compare Groups:**
   - Compare expression between controls and cases
   - Identify differentially expressed genes
   - Calculate fold changes and log2 fold changes
   - Visualize expression differences

**Key Analyses:**
- Sample counts by disease state
- Top expressed genes
- Genes upregulated/downregulated in cases vs. controls
- Expression statistics and distributions

### Best Practices for GEO Data Acquisition

1. **Dataset Selection:**
   - Choose datasets with appropriate sample sizes (typically n > 50)
   - Ensure phenotype data includes your variable of interest
   - Check that expression data is from relevant tissue types
   - Verify data quality and processing methods

2. **Data Quality:**
   - Check for batch effects (if multiple batches, consider batch correction)
   - Verify sample IDs match between expression and phenotype files
   - Check for missing values and outliers
   - Ensure expression data is properly normalized

3. **Platform Considerations:**
   - Different platforms (Affymetrix, Agilent, Illumina) have different probe designs
   - Probe-to-gene mappings may vary by platform
   - Some probes may not map to genes (control probes, intergenic regions)
   - Multiple probes per gene is common - aggregate appropriately

4. **Phenotype Extraction:**
   - Phenotype information may be in `Characteristics Ch1` or other columns
   - Parse structured data carefully (semicolon-separated key:value pairs)
   - Handle missing values appropriately
   - Create binary or categorical variables as needed for analysis

5. **File Management:**
   - Keep original downloaded files (Series Matrix, SOFT, MINiML) for reference
   - Save processed files (CSV) for analysis
   - Document which datasets were used and when they were downloaded
   - Note any data processing steps applied

### Troubleshooting

**Issue: Download fails with 404 error**
- Check that the GEO accession number is correct
- Verify the dataset is publicly available
- Try downloading manually from GEO website

**Issue: No gene mappings found**
- Check that the SOFT file contains platform annotation
- Verify the platform table section is being parsed correctly
- Some platforms may use different column names for gene information
- Check the SOFT file manually to see the actual structure

**Issue: Sample IDs don't match**
- Verify that sample IDs in expression data match phenotype data
- Check for different ID formats (with/without prefixes)
- Ensure both files are from the same GEO dataset

**Issue: Missing phenotype information**
- Check multiple columns in phenotype file (Description, Title, Characteristics)
- Some datasets may have limited phenotype data
- Consider downloading additional metadata files if needed

---

## Next Steps: Genotype Data Preparation

### Current Status ✅

**Completed:**
- ✅ GEO expression data downloaded and processed
- ✅ Phenotype data extracted and parsed
- ✅ Probe-to-gene mapping completed
- ✅ Chromosome positions extracted
- ✅ Data exploration and quality assessment

### Next Step: Download 1000 Genomes Genotype Data

**Objective:** Download VCF files from the 1000 Genomes Project to obtain genotype data for GWAS and eQTL analysis.

#### Step 1: Configure Genotype Download

Add to `config.yaml`:

```yaml
genotypes:
  source: "1000genomes"
  chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X"]
  # Or start with a subset for testing:
  # chromosomes: [2, 3, 4, 6, 8, 9, 10, 11]
  base_url: "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
  output_dir: "data/genotypes"
```

#### Step 2: Download VCF Files

**Implementation needed:**
- Create `src/download_genotypes.py` module
- Download VCF files for specified chromosomes
- Handle large file downloads with progress tracking
- Verify file integrity
- Support resume for interrupted downloads

**VCF File Format:**
- 1000 Genomes Phase 3 VCF files
- Format: `ALL.chr{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz`
- Each chromosome file is ~500MB - 2GB compressed

#### Step 3: Process VCF Files

**Tasks:**
- Parse VCF headers to extract sample IDs
- Filter variants (MAF, missingness, HWE)
- Convert to analysis-ready format
- Create sample ID mapping

### Important Considerations

#### Sample Matching Challenge

**Critical Issue:** GEO samples (GSM numbers) and 1000 Genomes samples are **different individuals**.

**Solutions:**

1. **For GWAS:**
   - Use 1000 Genomes samples with external phenotype data
   - Or use publicly available GWAS summary statistics
   - Or find datasets where the same individuals have both genotype and phenotype data

2. **For eQTL Analysis:**
   - Use resources like **GTEx** (Genotype-Tissue Expression) where samples have both genotype and expression
   - Or find datasets with matched genotype-expression pairs
   - Or use summary-level eQTL data from public databases

3. **Alternative Approach:**
   - Use GEO expression data for **differential expression analysis** (already done)
   - Use 1000 Genomes data for **population genetics** and **variant annotation**
   - Combine insights from both datasets at the interpretation level

#### Recommended Next Steps (Priority Order)

**Option A: Complete GWAS Pipeline (Recommended)**
1. Download 1000 Genomes genotype data
2. Use publicly available phenotype data or summary statistics
3. Perform GWAS on 1000 Genomes samples
4. Calculate PRS
5. Compare PRS with GEO expression patterns

**Option B: eQTL Analysis**
1. Download GTEx data (genotype + expression for same samples)
2. Or use summary-level eQTL data from public databases
3. Perform eQTL mapping
4. Integrate with GEO expression findings

**Option C: Integrative Analysis**
1. Use GEO data for expression analysis (disease-associated genes)
2. Use 1000 Genomes for variant annotation of those genes
3. Perform pathway analysis combining both datasets
4. Identify potential therapeutic targets

### Implementation Plan

**Immediate Next Step:** Download 1000 Genomes Genotype Data

**Files to Create:**
- `src/download_genotypes.py` - Download VCF files from 1000 Genomes
- `download_genotypes.py` - Main script to run downloads
- Update `config.yaml` with genotype settings

**After Genotype Download:**
- Quality control on genotype data
- Variant filtering
- Sample ID extraction
- Integration planning (based on chosen approach above)

---

## Data Exploration

See `notebooks/data_exploration.ipynb` for detailed exploration of downloaded GEO datasets.

---

## References

- GEO Database: https://www.ncbi.nlm.nih.gov/geo/
- SOFT Format Documentation: https://www.ncbi.nlm.nih.gov/geo/info/soft2.html
- MINiML Format Documentation: https://www.ncbi.nlm.nih.gov/geo/info/MINiML.html

