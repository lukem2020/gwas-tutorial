"""
GWAS Pipeline - Step-by-step implementation
"""

__version__ = "1.0.0"

from .config import Config, load_config
from .download_genotypes import GenotypeDownloader, download_genotypes

__all__ = ['Config', 'load_config', 'GenotypeDownloader', 'download_genotypes']

