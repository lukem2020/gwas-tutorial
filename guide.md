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

