#!/usr/bin/env python3
"""
Configuration loader for mbrain pipeline.

This module provides a centralized way to load configuration from config1.yaml.
All scripts should use this to get parameters instead of hardcoding values.

Usage:
    from config_loader import load_config, get_config
    
    # Load once at start (respects MBRAIN_PIPELINE env var)
    cfg = load_config()
    
    # Or explicitly specify pipeline
    cfg = load_config(pipeline='flip1')
    
    # Access parameters
    step_size = cfg['tractography']['step_size']
    max_dist = cfg['matlab_cleaning']['max_dist']
"""

import os
import yaml
from pathlib import Path
from typing import Dict, Any, Optional

# Global config cache
_config: Optional[Dict[str, Any]] = None
_current_pipeline: Optional[str] = None

def get_config_path(pipeline: Optional[str] = None) -> Path:
    """
    Get the path to config1.yaml.
    
    Returns:
        Path to config1.yaml
    """
    base_dir = Path(__file__).parent.parent
    return base_dir / "config1.yaml"

def load_config(config_path: Optional[Path] = None, pipeline: Optional[str] = None) -> Dict[str, Any]:
    """
    Load configuration from YAML file.
    
    Args:
        config_path: Optional explicit path to config file.
        pipeline: Ignored, kept for API compatibility.
        
    Returns:
        Dictionary with all configuration parameters.
    """
    global _config, _current_pipeline
    
    if config_path is None:
        config_path = get_config_path(pipeline)
    
    _current_pipeline = 'flip1'
    
    with open(config_path, 'r') as f:
        _config = yaml.safe_load(f)
    
    # Resolve {variable} placeholders in paths
    _resolve_path_variables(_config)
    
    # Allow environment variable override for derivatives directory
    # Used by --derivatives-dir CLI args across pipeline scripts
    deriv_override = os.environ.get('MBRAIN_DERIVATIVES_DIR')
    if deriv_override and 'paths' in _config:
        _config['paths']['derivatives_dir'] = deriv_override
        _config['paths']['derivatives_flip1'] = deriv_override
    
    return _config


def _resolve_path_variables(cfg: Dict[str, Any], max_passes: int = 5) -> None:
    """
    Resolve {variable} placeholders in the 'paths' section of the config.
    
    Variables like {repo_dir}, {data_base}, {matlab_path}, {rtp2_pipeline_path}
    are expanded using values from the same 'paths' dict.
    Also resolves {matlab_path} and {rtp2_pipeline_path} references in other sections.
    """
    if 'paths' not in cfg:
        return
    
    paths = cfg['paths']
    
    # Multi-pass resolution (handles chained references)
    for _ in range(max_passes):
        changed = False
        for key, val in paths.items():
            if isinstance(val, str) and '{' in val:
                new_val = val
                for ref_key, ref_val in paths.items():
                    if isinstance(ref_val, str) and '{' not in ref_val:
                        new_val = new_val.replace('{' + ref_key + '}', ref_val)
                if new_val != val:
                    paths[key] = new_val
                    changed = True
        if not changed:
            break
    
    # Also resolve path references in other sections (matlab_cleaning, atlas_coreg, file_tree, etc.)
    _resolve_nested(cfg, paths)


def _resolve_nested(obj: Any, paths: Dict[str, str]) -> None:
    """Recursively resolve {variable} placeholders in nested dicts/lists using paths dict."""
    if isinstance(obj, dict):
        for key, val in obj.items():
            if isinstance(val, str) and '{' in val:
                for ref_key, ref_val in paths.items():
                    if isinstance(ref_val, str):
                        val = val.replace('{' + ref_key + '}', ref_val)
                obj[key] = val
            elif isinstance(val, (dict, list)):
                _resolve_nested(val, paths)
    elif isinstance(obj, list):
        for i, item in enumerate(obj):
            if isinstance(item, str) and '{' in item:
                for ref_key, ref_val in paths.items():
                    if isinstance(ref_val, str):
                        item = item.replace('{' + ref_key + '}', ref_val)
                obj[i] = item
            elif isinstance(item, (dict, list)):
                _resolve_nested(item, paths)

def get_current_pipeline() -> str:
    """Get the currently loaded pipeline name ('flip1')."""
    global _current_pipeline
    if _current_pipeline is None:
        get_config()  # This will set _current_pipeline
    return _current_pipeline

def get_config() -> Dict[str, Any]:
    """
    Get the loaded configuration.
    
    Returns:
        Configuration dictionary. Loads from file if not already loaded.
    """
    global _config
    
    if _config is None:
        load_config()
    
    return _config

def get_tractography_params() -> Dict[str, Any]:
    """Get tractography (tckgen) parameters with backward-compatible aliases."""
    cfg = get_config()
    params = cfg['tractography'].copy()
    
    # Add backward-compatible aliases
    if 'step_size' in params and 'step' not in params:
        params['step'] = params['step_size']
    if 'min_length' in params and 'minlength' not in params:
        params['minlength'] = params['min_length']
    if 'max_length' in params and 'maxlength' not in params:
        params['maxlength'] = params['max_length']
    
    return params

def get_matlab_cleaning_params() -> Dict[str, Any]:
    """Get MATLAB cleaning parameters with backward-compatible aliases."""
    cfg = get_config()
    params = cfg['matlab_cleaning'].copy()
    
    # Add max_iter default if not present
    if 'max_iter' not in params:
        params['max_iter'] = 10
    
    return params

def get_ants_params() -> Dict[str, Any]:
    """Get ANTs registration parameters."""
    cfg = get_config()
    return cfg['ants']

def get_subject_groups() -> Dict[str, list]:
    """Get subject group definitions."""
    cfg = get_config()
    return cfg['groups']

def get_paths() -> Dict[str, str]:
    """Get all path configurations with backward-compatible aliases."""
    cfg = get_config()
    paths = cfg['paths'].copy()
    
    # Add backward-compatible aliases
    if 'repo_dir' in paths and 'base_dir' not in paths:
        paths['base_dir'] = paths['repo_dir']
    if 'allen_atlas_dir' in paths and 'atlas_dir' not in paths:
        paths['atlas_dir'] = paths['allen_atlas_dir']
    
    # Ensure atlas_rois_dir is set
    if 'atlas_rois_dir' not in paths:
        if 'allen_atlas_dir' in paths:
            paths['atlas_rois_dir'] = paths['allen_atlas_dir'] + '/rois'
        elif 'atlas_dir' in paths:
            paths['atlas_rois_dir'] = paths['atlas_dir'] + '/rois'
            
    if 'matlab_path' in paths and 'matlab_bin' not in paths:
        paths['matlab_bin'] = paths['matlab_path']
    if 'rtp2_pipeline_path' in paths and 'afq_path' not in paths:
        paths['afq_path'] = paths['rtp2_pipeline_path']
    
    return paths

def get_roi_definitions() -> Dict[str, Any]:
    """Get ROI definitions with Allen CCF IDs."""
    cfg = get_config()
    return cfg['rois']

def get_tract_definitions() -> list:
    """Get tract pair definitions."""
    cfg = get_config()
    return cfg['tracts']['pairs']

def get_manual_shift(subject_id: str) -> list:
    """
    Get manual atlas shift for a specific subject.
    
    Args:
        subject_id: Subject ID (e.g., 'sub-129-921-R1L1F')
        
    Returns:
        [x, y, z] shift in mm, or [0, 0, 0] if no shift defined.
    """
    cfg = get_config()
    shifts = cfg.get('atlas_coreg', {}).get('manual_shifts', {})
    return shifts.get(subject_id, [0.0, 0.0, 0.0])

def print_tractography_command(fod_path: str, output_path: str, 
                                seed_roi: str, include_roi: str,
                                mask_path: Optional[str] = None) -> str:
    """
    Generate tckgen command string using config parameters.
    
    Returns the command string for reference/logging.
    """
    params = get_tractography_params()
    
    cmd_parts = [
        "tckgen",
        fod_path,
        output_path,
        f"-algorithm {params['algorithm']}",
        f"-select {params['select']}",
        f"-step {params['step_size']}",
        f"-angle {params['angle']}",
        f"-minlength {params['min_length']}",
        f"-maxlength {params['max_length']}",
        f"-cutoff {params['cutoff']}",
        f"-seed_image {seed_roi}",
        f"-include {include_roi}",
    ]
    
    if mask_path and params.get('use_mask', True):
        cmd_parts.append(f"-mask {mask_path}")
    
    cmd_parts.append("-force")
    
    return " ".join(cmd_parts)

def get_matlab_cleaning_command(dwidir: str, tckfilename: str) -> str:
    """
    Generate MATLAB cleaning command string using config parameters.
    
    Returns the full matlab command for execution.
    """
    cfg = get_config()
    params = cfg['matlab_cleaning']
    paths = cfg['paths']
    
    matlab_cmd = (
        f"{paths['matlab_path']} -nodisplay -nosplash -r \""
        f"addpath(genpath('{paths['rtp2_pipeline_path']}')); "
        f"dwidir='{dwidir}'; "
        f"tckfilename='{tckfilename}'; "
        f"fg=fgRead(fullfile(dwidir,tckfilename)); "
        f"[fg_clean,~,SF]=AFQ_removeFiberOutliers(fg,"
        f"{params['max_dist']},{params['max_len']},{params['num_nodes']},"
        f"'{params['core_method']}',{params['verbose']},{params['max_iter']},{params['show_viz']}); "
        f"fg_SF=dtiNewFiberGroup([fg.name '_SF'],[0 0 155],[],[],SF.fibers); "
        f"AFQ_fgWrite(fg_clean,fullfile(dwidir,strrep(tckfilename,'.tck','_clean.tck')),'tck'); "
        f"AFQ_fgWrite(fg_SF,fullfile(dwidir,strrep(tckfilename,'.tck','_SF.tck')),'tck'); "
        f"fprintf('Done: %s -> %d fibers\\n',tckfilename,length(fg_clean.fibers)); "
        f"exit;\""
    )
    
    return matlab_cmd


if __name__ == "__main__":
    # Test loading config
    cfg = load_config()
    
    print("=== MBRAIN Configuration ===\n")
    
    print("Tractography parameters:")
    for k, v in get_tractography_params().items():
        print(f"  {k}: {v}")
    
    print("\nMATLAB cleaning parameters:")
    for k, v in get_matlab_cleaning_params().items():
        print(f"  {k}: {v}")
    
    print("\nSubject groups:")
    for group, subjects in get_subject_groups().items():
        print(f"  {group}: {subjects}")
    
    print("\nTract definitions:")
    if 'tracts' in cfg and 'pairs' in cfg['tracts']:
        for tract in cfg['tracts']['pairs']:
             print(f"  {tract.get('name', 'N/A')}: {tract.get('seed', 'N/A')} -> {tract.get('target', 'N/A')}")
    else:
        print("  No tract pairs defined in config.")
