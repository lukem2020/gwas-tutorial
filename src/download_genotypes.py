"""
Download genotype data from 1000 Genomes Project.
Step 2: Download VCF files for specified chromosomes.
"""

import os
import urllib.request
import urllib.error
from pathlib import Path
from typing import List, Optional
from tqdm import tqdm
import gzip
import shutil


class GenotypeDownloader:
    """Download 1000 Genomes genotype VCF files."""
    
    # 1000 Genomes Phase 3 base URL
    BASE_URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
    
    def __init__(self, output_dir: str, chromosomes: Optional[List[int]] = None):
        """
        Initialize genotype downloader.
        
        Parameters
        ----------
        output_dir : str
            Directory to save downloaded VCF files
        chromosomes : list of int, optional
            Chromosomes to download. If None, downloads all autosomes (1-22)
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        if chromosomes is None:
            chromosomes = list(range(1, 23))  # Autosomes
        self.chromosomes = chromosomes
    
    def get_vcf_filename(self, chromosome: int) -> str:
        """
        Get VCF filename for a chromosome.
        
        Parameters
        ----------
        chromosome : int
            Chromosome number
        
        Returns
        -------
        str
            VCF filename
        """
        return f"ALL.chr{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    
    def get_vcf_url(self, chromosome: int) -> str:
        """
        Get full URL for a chromosome VCF file.
        
        Parameters
        ----------
        chromosome : int
            Chromosome number
        
        Returns
        -------
        str
            Full FTP URL
        """
        filename = self.get_vcf_filename(chromosome)
        return f"{self.BASE_URL}/{filename}"
    
    def get_local_path(self, chromosome: int) -> Path:
        """
        Get local file path for a chromosome.
        
        Parameters
        ----------
        chromosome : int
            Chromosome number
        
        Returns
        -------
        Path
            Local file path
        """
        filename = self.get_vcf_filename(chromosome)
        return self.output_dir / filename
    
    def is_downloaded(self, chromosome: int) -> bool:
        """
        Check if chromosome VCF file is already downloaded.
        
        Parameters
        ----------
        chromosome : int
            Chromosome number
        
        Returns
        -------
        bool
            True if file exists and is not empty
        """
        path = self.get_local_path(chromosome)
        return path.exists() and path.stat().st_size > 0
    
    def download_chromosome(self, chromosome: int, force: bool = False) -> bool:
        """
        Download VCF file for a single chromosome.
        
        Parameters
        ----------
        chromosome : int
            Chromosome number
        force : bool
            If True, re-download even if file exists
        
        Returns
        -------
        bool
            True if download successful
        """
        if not force and self.is_downloaded(chromosome):
            print(f"  ✓ Chromosome {chromosome} already downloaded")
            return True
        
        url = self.get_vcf_url(chromosome)
        local_path = self.get_local_path(chromosome)
        temp_path = local_path.with_suffix('.tmp')
        
        try:
            print(f"  Downloading chromosome {chromosome}...")
            
            # Download with progress bar
            with urllib.request.urlopen(url) as response:
                total_size = int(response.headers.get('Content-Length', 0))
                
                with open(temp_path, 'wb') as f, tqdm(
                    total=total_size,
                    unit='B',
                    unit_scale=True,
                    unit_divisor=1024,
                    desc=f"    Chr{chromosome}",
                    miniters=1
                ) as pbar:
                    while True:
                        chunk = response.read(8192)
                        if not chunk:
                            break
                        f.write(chunk)
                        pbar.update(len(chunk))
            
            # Move temp file to final location
            temp_path.rename(local_path)
            
            file_size_mb = local_path.stat().st_size / (1024 * 1024)
            print(f"  ✓ Chromosome {chromosome} downloaded ({file_size_mb:.1f} MB)")
            return True
            
        except urllib.error.URLError as e:
            print(f"  ✗ Error downloading chromosome {chromosome}: {e}")
            if temp_path.exists():
                temp_path.unlink()
            return False
        except Exception as e:
            print(f"  ✗ Unexpected error for chromosome {chromosome}: {e}")
            if temp_path.exists():
                temp_path.unlink()
            return False
    
    def download_all(self, force: bool = False) -> dict:
        """
        Download VCF files for all specified chromosomes.
        
        Parameters
        ----------
        force : bool
            If True, re-download even if files exist
        
        Returns
        -------
        dict
            Summary with 'success', 'failed', and 'skipped' lists
        """
        print(f"Downloading 1000 Genomes Phase 3 VCF files...")
        print(f"Chromosomes: {self.chromosomes}")
        print(f"Output directory: {self.output_dir}")
        print()
        
        results = {
            'success': [],
            'failed': [],
            'skipped': []
        }
        
        for chromosome in self.chromosomes:
            if not force and self.is_downloaded(chromosome):
                results['skipped'].append(chromosome)
                continue
            
            if self.download_chromosome(chromosome, force=force):
                results['success'].append(chromosome)
            else:
                results['failed'].append(chromosome)
        
        return results
    
    def verify_downloads(self) -> dict:
        """
        Verify downloaded files exist and are not empty.
        
        Returns
        -------
        dict
            Verification results with 'valid' and 'invalid' lists
        """
        results = {
            'valid': [],
            'invalid': []
        }
        
        for chromosome in self.chromosomes:
            path = self.get_local_path(chromosome)
            if path.exists() and path.stat().st_size > 0:
                results['valid'].append(chromosome)
            else:
                results['invalid'].append(chromosome)
        
        return results


def download_genotypes(config, force: bool = False) -> dict:
    """
    Download 1000 Genomes genotype data using configuration.
    
    Parameters
    ----------
    config
        Configuration object
    force : bool
        If True, re-download even if files exist
    
    Returns
    -------
    dict
        Download results summary
    """
    output_dir = config.get('paths.genotypes_dir')
    chromosomes = config.get('genotypes.chromosomes')
    
    downloader = GenotypeDownloader(output_dir, chromosomes)
    results = downloader.download_all(force=force)
    
    return results

