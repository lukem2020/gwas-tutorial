"""
Configuration management for GWAS pipeline.
Step 1: Load and validate project configuration.
"""

import yaml
from pathlib import Path
from typing import Dict, Any, Optional
import os


class Config:
    """Configuration manager for GWAS pipeline."""
    
    def __init__(self, config_path: str = "config.yaml"):
        """
        Initialize configuration from YAML file.
        
        Parameters
        ----------
        config_path : str
            Path to configuration YAML file
        """
        self.config_path = Path(config_path)
        self.config = self._load_config()
        self._validate_config()
        self._create_directories()
    
    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from YAML file."""
        if not self.config_path.exists():
            raise FileNotFoundError(
                f"Configuration file not found: {self.config_path}\n"
                f"Please create config.yaml in the project root."
            )
        
        with open(self.config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        return config
    
    def _validate_config(self) -> None:
        """Validate configuration structure and values."""
        required_sections = [
            'project', 'paths', 'genotypes', 'phenotypes', 
            'qc', 'pca', 'gwas'
        ]
        
        for section in required_sections:
            if section not in self.config:
                raise ValueError(f"Missing required configuration section: {section}")
        
        # Validate paths
        if 'data_dir' not in self.config['paths']:
            raise ValueError("Missing 'data_dir' in paths configuration")
        
        # Validate QC thresholds
        qc = self.config['qc']
        if not (0 < qc['sample_missingness_threshold'] <= 1):
            raise ValueError("sample_missingness_threshold must be between 0 and 1")
        if not (0 < qc['maf_threshold'] < 0.5):
            raise ValueError("maf_threshold must be between 0 and 0.5")
    
    def _create_directories(self) -> None:
        """Create necessary directories if they don't exist."""
        paths = self.config['paths']
        directories = [
            paths['data_dir'],
            paths['genotypes_dir'],
            paths['phenotypes_dir'],
            paths['expression_dir'],
            paths['results_dir'],
            paths['plots_dir']
        ]
        
        for directory in directories:
            Path(directory).mkdir(parents=True, exist_ok=True)
    
    def get(self, key: str, default: Any = None) -> Any:
        """
        Get configuration value using dot notation.
        
        Parameters
        ----------
        key : str
            Configuration key in dot notation (e.g., 'qc.maf_threshold')
        default : Any
            Default value if key not found
        
        Returns
        -------
        Any
            Configuration value
        """
        keys = key.split('.')
        value = self.config
        
        for k in keys:
            if isinstance(value, dict) and k in value:
                value = value[k]
            else:
                return default
        
        return value
    
    def __getitem__(self, key: str) -> Any:
        """Allow dictionary-like access to configuration."""
        return self.config[key]
    
    def __contains__(self, key: str) -> bool:
        """Check if configuration key exists."""
        return key in self.config
    
    def save(self, output_path: Optional[str] = None) -> None:
        """
        Save current configuration to file.
        
        Parameters
        ----------
        output_path : str, optional
            Path to save configuration. If None, saves to original path.
        """
        path = Path(output_path) if output_path else self.config_path
        
        with open(path, 'w') as f:
            yaml.dump(self.config, f, default_flow_style=False, sort_keys=False)
    
    def update(self, key: str, value: Any) -> None:
        """
        Update configuration value using dot notation.
        
        Parameters
        ----------
        key : str
            Configuration key in dot notation
        value : Any
            New value
        """
        keys = key.split('.')
        config = self.config
        
        for k in keys[:-1]:
            if k not in config:
                config[k] = {}
            config = config[k]
        
        config[keys[-1]] = value


def load_config(config_path: str = "config.yaml") -> Config:
    """
    Convenience function to load configuration.
    
    Parameters
    ----------
    config_path : str
        Path to configuration file
    
    Returns
    -------
    Config
        Configuration object
    """
    return Config(config_path)


if __name__ == "__main__":
    # Test configuration loading
    try:
        config = load_config()
        print("✓ Configuration loaded successfully!")
        print(f"  Project: {config.get('project.name')}")
        print(f"  Data directory: {config.get('paths.data_dir')}")
        print(f"  Chromosomes: {config.get('genotypes.chromosomes')}")
        print(f"  QC MAF threshold: {config.get('qc.maf_threshold')}")
    except Exception as e:
        print(f"✗ Error loading configuration: {e}")

