"""
main.py is the entry point for the pipeline.
It downloads the GEO datasets, aggregates the files, and cleans the data.

Modules used:
- download_geo_data.py: downloads the GEO datasets
- aggregate_files.py: aggregates the files
- clean_data.py: cleans the data

Configuration:
- config.yaml: configuration for the pipeline

Data directories:
- data/expression: directory for the aggregated files
- data/expression/clean: directory for the cleaned data

Output:
- data/expression/clean: directory for the cleaned data
"""
from src.download_geo_data import download_all_geo_datasets
from src.aggregate_files import process_all_datasets
from src.clean_data import clean_all_datasets


if __name__ == "__main__":
    download_all_geo_datasets()
    process_all_datasets()
    clean_all_datasets()