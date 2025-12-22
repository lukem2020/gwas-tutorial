> **Step-by-step Python GWAS+PRS+eQTL Tutorial (1000 Genomes Real Data)**

This guide covers a *full pipeline* for GWAS, polygenic risk scoring, and eQTL analysis using *real genotype data* from the 1000 Genomes Project. All steps can be run in pure Python (cross-platform, no external bioinformatics tools required).

---

## Pipeline Overview

| Step | Task                      | Time       | Output               |
|------|---------------------------|------------|----------------------|
| 1    | Configure project         | < 1 sec    | config.yaml          |
| 2    | Download genotype data    | 10-60 min  | VCF files (chr 2,3,4,6,8,9,10,11) |
| 3    | Prepare phenotype data   | 1-5 min    | CSV file (optional)  |
| 4    | Load genotypes            | 1-30 min   | Genotype matrix      |
| 5    | Load phenotype            | < 1 sec    | Cases/controls       |
| 6    | Run GWAS                  | 3-30 min   | P-values, Manhattan plot |
| 7    | Calculate PRS             | 2-5 min    | AUC scores, quantiles |
| 8    | Download expression data  | 5-15 min   | GEO expression matrix (optional) |
| 9    | Run eQTL analysis         | 5-60 min   | SNP-gene associations (optional) |
| 10   | Visualize results         | 5-10 sec   | Plots, summary stats |
| 11   | Save results              | 1-2 sec    | CSV files, plots     |
| 12   | Review & QC               | 10-30 min  | Quality checks       |

---

## GWAS Best Practices

### 1. Quality Control (QC) - Critical First Steps

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

### 2. Statistical Considerations

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

### 3. Data Preprocessing

#### Genotype Data
- **Imputation**: Impute missing genotypes using reference panels (1000 Genomes, HRC, TOPMed)
- **Phasing**: Phase haplotypes for better imputation accuracy
- **Linkage Disequilibrium (LD)**: Calculate LD structure for regional analysis
- **Allele coding**: Ensure consistent allele coding (reference vs. alternate)

#### Phenotype Data
- **Trait transformation**: Consider log/normalization for skewed continuous traits
- **Outlier removal**: Remove extreme outliers that may be errors
- **Covariate adjustment**: Pre-adjust phenotypes for known confounders if needed

### 4. Population Stratification Control

#### Principal Component Analysis (PCA)
- Calculate PCs on LD-pruned variants (r² < 0.2)
- Include top PCs (typically 3-10) as covariates
- Visualize PC plots to identify population clusters
- Consider stratified analysis if distinct populations present

#### Alternative Methods
- **Structured association**: Use STRUCTURE/ADMIXTURE for ancestry inference
- **Mixed models**: Use linear mixed models (LMM) to account for relatedness
- **Stratified analysis**: Analyze populations separately if highly heterogeneous

### 5. Reproducibility & Documentation

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

### 6. Visualization Best Practices

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

### 7. Interpretation Guidelines

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

### 8. Common Pitfalls to Avoid

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

### 9. Post-GWAS Analysis

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

### 10. Ethical Considerations

- **Data privacy**: Ensure proper consent and data protection
- **Population representation**: Acknowledge limitations in diverse populations
- **Clinical interpretation**: Be cautious about clinical translation
- **Reporting**: Avoid deterministic language about genetic risk

---

## Recommended Workflow Checklist

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

