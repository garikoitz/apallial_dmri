"""
MBRAIN Functions Module

This module contains all Python scripts for the mouse brain tractography pipeline.
All scripts use config_loader to read parameters from config1.yaml.

Usage:
    from functions.config_loader import load_config
    from functions import atlas_coregistration
"""

from .config_loader import load_config, get_config, get_paths, get_subject_groups
