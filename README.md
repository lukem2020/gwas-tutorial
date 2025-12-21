# gwas-tutorial

GWAS tutorial project for working with genetic data using Python-only packages (works on Windows, Mac, and Linux).

## Setup

1. Create and activate a virtual environment:
```bash
python -m venv .venv
.venv\Scripts\activate  # Windows
# or
source .venv/bin/activate  # Mac/Linux
```

2. Install required packages:
```bash
pip install -r requirements.txt
```

## Working with VCF Files

This project uses **scikit-allel**, a pure Python package that works on Windows without requiring external tools like tabix or bcftools.

### Reading VCF Files

```python
from src.utils.vcf_indexer import read_vcf

# Read entire VCF file
callset = read_vcf('data/1000genomes/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz')

# Read specific region
callset = read_vcf('data/file.vcf.gz', region='22:1000000-2000000')
```

### Using the Command Line

```bash
# Read VCF file metadata
python -m src.utils.vcf_indexer data/1000genomes/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

# Read specific region
python -m src.utils.vcf_indexer data/1000genomes/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz 22:1000000-2000000
```

## Python Packages Used

- **scikit-allel**: VCF file reading and genetic data manipulation (pure Python, works on Windows)
- **pandas**: Data manipulation
- **numpy**: Numerical computing
- **scipy**: Statistical functions
- **statsmodels**: Statistical modeling
- **scikit-learn**: Machine learning (for PRS calculations)
- **matplotlib**: Plotting and visualization

All packages are pure Python or have Windows wheels available - no compilation needed!