# Windows Setup for GWAS Tutorial

## The Problem

Many bioinformatics tools like `pysam`, `tabix`, and `bcftools` are not available for Windows via conda. They require Unix-like environments.

## Solutions

### Option 1: Use WSL (Windows Subsystem for Linux) - **Recommended**

1. **Install WSL** (if not already installed):
   ```powershell
   wsl --install
   ```

2. **In WSL, set up your environment:**
   ```bash
   # Navigate to your project (mount Windows drive)
   cd /mnt/c/Users/User/Desktop/Personal/github/projects/gwas/gwas-tutorial
   
   # Create conda environment
   conda create -n gwas-tutorial python=3.11 -y
   conda activate gwas-tutorial
   conda install -c bioconda tabix bcftools pysam -y
   ```

3. **Index your VCF file:**
   ```bash
   tabix -p vcf data/1000genomes/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
   ```

### Option 2: Download Pre-built Windows Binaries

1. **Download HTSlib for Windows:**
   - Visit: https://github.com/samtools/htslib/releases
   - Download the Windows binary (if available) or build from source
   - Extract and add to PATH

2. **Or use pre-built binaries from:**
   - https://github.com/samtools/samtools/releases (includes tabix)

### Option 3: Use Python-only Alternatives

For VCF indexing, you could use Python libraries that don't require external tools, though they may be slower:
- `cyvcf2` (if it installs)
- Pure Python implementations

### Option 4: Use Docker

Run a Linux container with all tools pre-installed:
```bash
docker run -it -v /path/to/data:/data biocontainers/samtools
```

## Current Status

- ✅ Python environment works
- ❌ `pysam` - Not available for Windows via conda
- ❌ `tabix` - Not available for Windows via conda  
- ❌ `bcftools` - Not available for Windows via conda

**Recommendation:** Use WSL for the best experience with bioinformatics tools on Windows.


