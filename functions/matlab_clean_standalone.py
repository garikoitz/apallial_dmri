#!/usr/bin/env python3
"""
Standalone MATLAB tract cleaning script.
Designed to run without terminal interaction.
"""

import subprocess
import sys
from pathlib import Path
from datetime import datetime

from config_loader import load_config

def get_tracts_needing_cleaning():
    cfg = load_config()
    deriv_dir = Path(cfg['paths']['derivatives_dir'])
    
    tract_patterns = [
        "tracks_Isocortex_Striatum_*.tck",
        "tracks_Isocortex_Thalamus_*.tck",
        "tracks_Graft_Striatum_*.tck",
        "tracks_Graft_Thalamus_*.tck",
    ]
    
    tracts = []
    for pattern in tract_patterns:
        for tck_file in deriv_dir.rglob(pattern):
            if any(x in str(tck_file) for x in ['_cleaned', '_sift', '_top', '_SF', '_warped', '_both_rois', '_clean.']):
                continue
            if '/old/' in str(tck_file) or '/template/' in str(tck_file) or '/tracts/' in str(tck_file):
                continue
            
            cleaned_file = tck_file.parent / tck_file.name.replace('.tck', '_cleaned.tck')
            if not cleaned_file.exists():
                tracts.append((tck_file, cleaned_file))
    
    return sorted(tracts, key=lambda x: str(x[0]))

def clean_tract(input_file, output_file, cfg, log_file):
    matlab_cfg = cfg.get('matlab_cleaning', {})
    matlab_bin = matlab_cfg.get('matlab_bin', '/Applications/MATLAB_R2025b.app/bin/matlab')
    afq_path = matlab_cfg.get('afq_path', '/Users/glerma/soft/gitlab/rtp2-pipeline')
    max_dist = matlab_cfg.get('max_dist', 3)
    max_len = matlab_cfg.get('max_len', 3)
    num_nodes = matlab_cfg.get('num_nodes', 100)
    core_method = matlab_cfg.get('core_method', 'mean')
    
    dwi_dir = input_file.parent
    sf_file = dwi_dir / input_file.name.replace('.tck', '_SF.tck')
    
    # Build MATLAB command - filter short fibers first, then run AFQ
    matlab_cmd = (
        f"addpath(genpath('{afq_path}')); "
        f"dwidir = '{dwi_dir}'; "
        f"tckfilename = '{input_file.name}'; "
        f"try, "
        f"fg = fgRead(fullfile(dwidir, tckfilename)); "
        f"nOrig = length(fg.fibers); "
        f"lens = cellfun(@(f) size(f,2), fg.fibers); "
        f"fg.fibers = fg.fibers(lens > 3); "  # Filter fibers with <=3 points
        f"fprintf('Input: %s with %d fibers (filtered %d short)\\n', tckfilename, length(fg.fibers), nOrig - length(fg.fibers)); "
        f"[fg_clean, ~, SF] = AFQ_removeFiberOutliers(fg, {max_dist}, {max_len}, {num_nodes}, '{core_method}', 0, 10, 0); "
        f"fg_SF = dtiNewFiberGroup([fg.name '_SF'], [0 0 155], [], [], SF.fibers); "
        f"AFQ_fgWrite(fg_clean, fullfile(dwidir, '{output_file.name}'), 'tck'); "
        f"AFQ_fgWrite(fg_SF, fullfile(dwidir, '{sf_file.name}'), 'tck'); "
        f"fprintf('Output: %s with %d fibers (removed %d)\\n', '{output_file.name}', length(fg_clean.fibers), nOrig - length(fg_clean.fibers)); "
        f"catch ME, "
        f"fprintf('ERROR: %s\\n', ME.message); "
        f"end; "
        f"exit;"
    )
    
    cmd = [matlab_bin, "-nodisplay", "-nosplash", "-batch", matlab_cmd]
    
    with open(log_file, 'a') as f:
        f.write(f"\n[{datetime.now()}] Processing: {input_file.name}\n")
        f.flush()
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            
            if result.returncode != 0 or not output_file.exists():
                error_msg = result.stderr[:500] if result.stderr else result.stdout[:500] if result.stdout else 'Unknown'
                f.write(f"  ERROR: {error_msg}\n")
                return False
            
            f.write(f"  SUCCESS: {output_file.name} created\n")
            return True
            
        except subprocess.TimeoutExpired:
            f.write(f"  ERROR: Timeout after 600s (10 min)\n")
            return False
        except Exception as e:
            f.write(f"  ERROR: {str(e)}\n")
            return False

def main():
    log_file = f"/tmp/matlab_clean_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    cfg = load_config()
    
    with open(log_file, 'w') as f:
        f.write(f"MATLAB Tract Cleaning - Started {datetime.now()}\n")
        f.write("="*60 + "\n")
    
    tracts = get_tracts_needing_cleaning()
    
    with open(log_file, 'a') as f:
        f.write(f"Tracts to clean: {len(tracts)}\n\n")
    
    success = 0
    for i, (input_file, output_file) in enumerate(tracts, 1):
        print(f"[{i}/{len(tracts)}] {input_file.name}...", flush=True)
        if clean_tract(input_file, output_file, cfg, log_file):
            success += 1
            print(f"  ✓ Done")
        else:
            print(f"  ✗ Failed")
    
    with open(log_file, 'a') as f:
        f.write(f"\nCompleted: {success}/{len(tracts)} tracts cleaned\n")
        f.write(f"Finished: {datetime.now()}\n")
    
    print(f"\nDone: {success}/{len(tracts)} cleaned. Log: {log_file}")
    return 0 if success == len(tracts) else 1

if __name__ == '__main__':
    sys.exit(main())
