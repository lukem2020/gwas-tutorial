# gwas-tutorial

GWAS tutorial project for working with genetic data.

## Setup

### Using Conda (Recommended for Windows)

1. Create and activate the conda environment:
```bash
conda env create -f environment.yml
conda activate gwas-tutorial
```

2. If you already have a conda environment, install the packages:
```bash
conda install -c conda-forge -c bioconda pysam tabix bcftools
```

### Indexing VCF Files

To index a VCF file, you can use either:

**Option 1: Using tabix command (after conda install)**
```bash
tabix -p vcf data/1000genomes/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
```

**Option 2: Using Python script**
```bash
python -m src.utils.vcf_indexer data/1000genomes/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
```