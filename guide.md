> **Step-by-step Python GWAS+PRS+eQTL Tutorial (1000 Genomes & GEO Data)**

This guide covers a *full pipeline* for GWAS, polygenic risk scoring, and eQTL analysis using *real genotype data* from the **1000 Genomes Project** and *expression data* from **GEO (Gene Expression Omnibus)**. All steps can be run in pure Python (cross-platform, no external bioinformatics tools required).

---

## Pipeline Overview

### Pre-Analysis: Data Preparation & Quality Control

| Step | Task                      | Time       | Output               |
|------|---------------------------|------------|----------------------|
| 1    | Configure project         | < 1 sec    | config.yaml          |
| 2    | Download genotype data (1000 Genomes) | 10-60 min  | VCF files (chr 2,3,4,6,8,9,10,11) |
| 3    | Prepare phenotype data   | 1-5 min    | CSV file             |
| 4    | Download expression data (GEO) | 5-15 min   | GEO expression matrix (optional, for eQTL) |
| 5    | Load genotypes            | 1-30 min   | Genotype matrix      |
| 6    | Load phenotype            | < 1 sec    | Cases/controls       |
| 7    | Quality Control (QC)      | 10-30 min  | Filtered samples/variants, QC metrics |
| 8    | Calculate PCA             | 2-10 min   | Principal components, population structure plots |

### Post-Analysis: Association Testing & Interpretation

| Step | Task                      | Time       | Output               |
|------|---------------------------|------------|----------------------|
| 9    | Run GWAS                  | 3-30 min   | P-values, effect sizes, Manhattan plot |
| 10   | Calculate PRS             | 2-5 min    | PRS scores, AUC, quantiles |
| 11   | Run eQTL analysis         | 5-60 min   | SNP-gene associations (optional) |
| 12   | Visualize results         | 5-10 sec   | Manhattan, QQ plots, regional plots |
| 13   | Save results              | 1-2 sec    | CSV files, plots, summary stats |
| 14   | Review & interpret        | 10-30 min  | Annotated results, validation |

---

## GWAS Best Practices

### 1. Data Acquisition

> **Note**: This tutorial uses **1000 Genomes Project** for genotype data and **GEO** for expression data. These are complementary public datasets that enable GWAS and eQTL analysis workflows.

#### Genotype Data Sources
- **Primary source**: **1000 Genomes Project** (https://www.internationalgenome.org/)
  - Phase 3 release: ~2,500 individuals from 26 populations
  - VCF format available for download
  - Includes population labels and sample metadata
- **Format considerations**: VCF format (can be converted to PLINK bed/bim/fam, HDF5, or Zarr if needed)
- **Data versioning**: Document which 1000 Genomes phase/version you're using (e.g., Phase 3, GRCh37/hg19)
- **Download verification**: Check file integrity (MD5/SHA checksums provided by 1000 Genomes)
- **Storage**: Plan for large file sizes (genome-wide VCF files can be 100s of GB per chromosome)
- **Access**: Publicly available, no restrictions for research use

#### Phenotype Data
- **Clinical data**: Ensure proper consent and IRB approval
- **Data formats**: CSV, TSV, or database formats
- **Data dictionary**: Document all variables, coding schemes, and units
- **Missing data**: Document missing data patterns and reasons
- **Trait definition**: Clearly define case/control criteria or continuous trait measurements

#### Expression Data (for eQTL Analysis)
- **Primary source**: **GEO (Gene Expression Omnibus)** (https://www.ncbi.nlm.nih.gov/geo/)
  - Largest repository of gene expression data
  - Contains microarray and RNA-seq datasets
  - Searchable by organism, tissue, disease, platform
  - Access via GEO accession numbers (e.g., GSE12345)
  - Download processed data matrices or raw files
- **GEO data formats**: 
  - **Series Matrix files**: Pre-processed expression matrices (recommended for quick start)
  - **Raw data**: CEL files (Affymetrix), IDAT files (Illumina), FASTQ (RNA-seq)
  - **Processed data**: Gene expression matrices (TPM, FPKM, RPKM, counts, normalized intensities)
  - **Standard formats**: GCT, SOFT format files
- **Data processing considerations**:
  - **Normalization**: RPKM, FPKM, TPM for RNA-seq; RMA, quantile normalization for microarrays
  - **Batch effects**: Document and account for batch/plate effects
  - **Quality metrics**: RIN scores (RNA-seq), detection p-values (microarrays)
  - **Gene annotation**: Use consistent gene ID systems (Ensembl, RefSeq, Gene Symbol)
- **Sample matching for eQTL analysis**: 
  - **Direct matching**: Ideally use GEO datasets where expression samples correspond to 1000 Genomes individuals (check GEO sample metadata for sample IDs)
  - **Population-matched analysis**: If direct matching isn't available, use population-matched expression data (e.g., European samples from GEO with European samples from 1000 Genomes)
  - **Cross-reference workflow**: Match by population labels and perform population-stratified eQTL analysis
  - **Documentation**: Always document sample overlap, matching strategy, and any assumptions made
- **Tissue/cell type**: Document tissue type, cell line, or experimental condition
- **Metadata**: Document experimental design, protocols, and processing pipelines
- **Storage**: Expression matrices can be large (especially RNA-seq with many samples/genes)

#### Reference Data
- **Reference panels**: **1000 Genomes** (primary reference for this tutorial), HRC, TOPMed for imputation
- **Annotation databases**: dbSNP, Ensembl, RefSeq for variant annotation
- **LD reference**: Use 1000 Genomes LD structure for population-specific analysis
- **Functional annotation**: ENCODE, Roadmap Epigenomics for functional elements
- **GEO metadata**: Use GEO sample and platform annotations for expression data interpretation

#### Data Organization
- **File naming**: Use consistent, descriptive naming conventions
- **Directory structure**: Organize by data type (genotypes, phenotypes, results)
- **Metadata**: Create README files documenting data sources and processing
- **Backup**: Ensure proper backup and version control for data files

### 2. Quality Control (QC) - Critical First Steps

#### Sample-Level QC
- **Missingness**: Remove samples with >5-10% missing genotype data
- **Sex check**: Verify reported sex matches genetic sex (X chromosome heterozygosity)
- **Relatedness**: Remove or account for related individuals (pi-hat > 0.2)
- **Population stratification**: Use PCA to identify and control for population structure
- **Heterozygosity**: Remove outliers (excessive heterozygosity may indicate contamination)

#### Variant-Level QC
- **Missingness**: Remove variants with >5-10% missing data
- **Minor Allele Frequency (MAF)**: Filter rare variants (typically MAF < 0.01 or 0.05)
- **Hardy-Weinberg Equilibrium (HWE)**: Remove variants deviating from HWE (p < 1e-6 in controls)
- **Call rate**: Ensure high-quality genotyping (call rate > 95-99%)
- **Duplicate variants**: Remove duplicate or ambiguous variants

### 3. Statistical Considerations

#### Association Testing
- **Model selection**: Use appropriate regression model (logistic for binary traits, linear for continuous)
- **Covariates**: Include age, sex, and principal components (PCs) as covariates
- **Population stratification**: Include top 3-10 principal components to control for ancestry
- **Genomic control**: Check λ (genomic inflation factor) - should be close to 1.0
  - λ > 1.1 suggests population stratification or other confounders
  - λ < 1.0 suggests over-correction

#### Multiple Testing Correction
- **Bonferroni correction**: Conservative threshold (p < 5×10⁻⁸ for genome-wide significance)
- **False Discovery Rate (FDR)**: Less conservative alternative (Benjamini-Hochberg)
- **Permutation testing**: Gold standard for empirical p-values
- **Regional significance**: Consider significance within LD blocks

### 4. Data Preprocessing

#### Genotype Data
- **Imputation**: Impute missing genotypes using **1000 Genomes** as reference panel (HRC, TOPMed are alternatives)
- **Phasing**: Phase haplotypes for better imputation accuracy (1000 Genomes provides phased data)
- **Linkage Disequilibrium (LD)**: Calculate LD structure from 1000 Genomes data for regional analysis
- **Allele coding**: Ensure consistent allele coding (reference vs. alternate) - 1000 Genomes uses GRCh37/hg19 reference

#### Phenotype Data
- **Trait transformation**: Consider log/normalization for skewed continuous traits
- **Outlier removal**: Remove extreme outliers that may be errors
- **Covariate adjustment**: Pre-adjust phenotypes for known confounders if needed

### 5. Population Stratification Control

#### Principal Component Analysis (PCA)
- Calculate PCs on LD-pruned variants (r² < 0.2)
- Include top PCs (typically 3-10) as covariates
- Visualize PC plots to identify population clusters
- Consider stratified analysis if distinct populations present

#### Alternative Methods
- **Structured association**: Use STRUCTURE/ADMIXTURE for ancestry inference
- **Mixed models**: Use linear mixed models (LMM) to account for relatedness
- **Stratified analysis**: Analyze populations separately if highly heterogeneous

### 6. Reproducibility & Documentation

#### Code & Data Management
- **Version control**: Use Git for code versioning
- **Environment**: Document Python/R versions and package versions (requirements.txt)
- **Random seeds**: Set random seeds for reproducibility
- **Configuration files**: Use YAML/JSON configs for parameters
- **Data provenance**: Document data sources, versions, and processing steps

#### Results Documentation
- **QC metrics**: Report all QC filters and thresholds applied
- **Sample sizes**: Report final sample sizes after QC
- **Effect sizes**: Report odds ratios/beta coefficients with confidence intervals
- **Manhattan plots**: Include genome-wide association plots
- **QQ plots**: Include quantile-quantile plots to assess inflation

### 7. Visualization Best Practices

#### Essential Plots
- **Manhattan plot**: Show -log10(p-values) across chromosomes
- **QQ plot**: Assess genomic inflation and identify systematic biases
- **Effect size plots**: Visualize odds ratios/beta coefficients
- **Regional plots**: Zoom into significant regions with LD information
- **PCA plots**: Show population structure

#### Plot Guidelines
- Use genome-wide significance threshold (p = 5×10⁻⁸) as horizontal line
- Use suggestive threshold (p = 1×10⁻⁵) as secondary line
- Color-code by chromosome for Manhattan plots
- Include sample sizes and QC metrics in plot titles/captions

### 8. Interpretation Guidelines

#### Significance Assessment
- **Genome-wide significant**: p < 5×10⁻⁸ (gold standard)
- **Suggestive**: p < 1×10⁻⁵ (may warrant replication)
- **Replication**: Significant findings should be replicated in independent cohorts
- **Effect sizes**: Consider clinical/biological significance, not just statistical

#### Biological Interpretation
- **Gene annotation**: Identify nearest genes and functional elements
- **Functional annotation**: Check if variants are in coding regions, regulatory elements
- **Pathway analysis**: Test for enrichment in biological pathways
- **Literature review**: Compare with previous GWAS findings

### 9. Common Pitfalls to Avoid

#### Data Issues
- ❌ Not checking for batch effects or plate effects
- ❌ Ignoring population stratification
- ❌ Using raw p-values without multiple testing correction
- ❌ Not filtering low-quality variants/samples
- ❌ Mixing populations without proper control

#### Analysis Issues
- ❌ Overfitting models with too many covariates
- ❌ Not accounting for relatedness in the sample
- ❌ Ignoring genomic inflation (λ)
- ❌ Not validating findings in independent data
- ❌ Interpreting suggestive hits as definitive

#### Reporting Issues
- ❌ Not reporting effect sizes, only p-values
- ❌ Not documenting QC steps
- ❌ Over-interpreting marginal associations
- ❌ Not considering LD structure in interpretation

### 10. Post-GWAS Analysis

#### Fine-Mapping
- Identify causal variants within significant regions
- Use statistical fine-mapping (e.g., FINEMAP, SuSiE)
- Consider functional annotation to prioritize variants

#### Functional Follow-up
- **eQTL analysis**: Test if variants affect gene expression
- **Colocalization**: Test if GWAS and eQTL signals colocalize
- **Mendelian randomization**: Test causal relationships
- **Functional experiments**: Validate findings experimentally

#### Polygenic Risk Scores (PRS)
- Use independent discovery and target datasets
- Optimize p-value thresholds
- Validate PRS in independent cohorts
- Report predictive performance (AUC, R²)

### 11. Ethical Considerations

- **Data privacy**: Ensure proper consent and data protection
- **Population representation**: Acknowledge limitations in diverse populations
- **Clinical interpretation**: Be cautious about clinical translation
- **Reporting**: Avoid deterministic language about genetic risk

---

## Recommended Workflow Checklist

- [ ] Acquire and verify 1000 Genomes genotype data (check integrity, document phase/version)
- [ ] Acquire and prepare phenotype data (document variables, handle missing data)
- [ ] Acquire expression data from GEO (if doing eQTL analysis: search GEO, download Series Matrix or processed data)
- [ ] Organize data with clear directory structure and metadata
- [ ] Perform sample-level QC (missingness, relatedness, sex check)
- [ ] Perform variant-level QC (MAF, HWE, missingness)
- [ ] Calculate and visualize PCA
- [ ] Include appropriate covariates (age, sex, PCs)
- [ ] Run association tests with proper model
- [ ] Check genomic inflation (λ)
- [ ] Apply multiple testing correction
- [ ] Generate Manhattan and QQ plots
- [ ] Document all QC steps and thresholds
- [ ] Validate significant findings (replication if possible)
- [ ] Annotate significant variants
- [ ] Report effect sizes with confidence intervals

