#!/usr/bin/env python3
"""
Collect bilateral passthrough tract data for report generation.
"""

import os
import subprocess
import json
import sys

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DERIV = "/Users/glerma/Library/CloudStorage/GoogleDrive-garikoitz@gmail.com/Shared drives/DTI Collab/mbrain_data/data/derivatives/mrtrix_preproc_flip1"

subjects = {
    'CTRL': ['sub-129-921-R1F', 'sub-129-921-R1L1F', 'sub-129-1014-L1M', 'sub-138-1010-L1F'],
    'AP': ['sub-596-C-OO-NmM', 'sub-769-C-O5-R1M', 'sub-769-C-OO-NmM', 'sub-875-C-OO-NM'],
    'APTR': ['sub-AP72-862-901-NmF', 'sub-AP72-862-901-R1M', 'sub-Apallial17NmF', 'sub-Apallial20R1F']
}

ctrl_segments = ['anterior', 'anterior_central', 'central', 'posterior_central', 'posterior']
graft_segments = ['anterior_dorsal', 'central_dorsal', 'posterior_dorsal']
targets = {
    'CTRL': ['thalamus', 'striatum', 'cerebellum', 'olfactory_cortex'],
    'AP': ['Thalamus', 'Striatum', 'cerebellum', 'Palecortex'],
    'APTR': ['Thalamus', 'Striatum', 'cerebellum', 'Palecortex']
}

def get_session(subj_path):
    """Get the session folder for a subject."""
    try:
        sessions = [d for d in os.listdir(subj_path) if d.startswith('ses-')]
        return sessions[0] if sessions else None
    except:
        return None

def get_tck_count(tck_path):
    """Get streamline count from a .tck file."""
    try:
        result = subprocess.run(['tckinfo', tck_path], capture_output=True, text=True)
        for line in result.stdout.split('\n'):
            if 'count:' in line and 'total_count' not in line:
                return int(line.split(':')[1].strip())
    except:
        pass
    return 0

def collect_hemisphere_data(hemisphere='right'):
    """Collect tract data for one hemisphere."""
    hemi = hemisphere
    data = {'CTRL': {}, 'AP': {}, 'APTR': {}}
    
    for group, subj_list in subjects.items():
        for subj in subj_list:
            subj_path = os.path.join(DERIV, subj)
            ses = get_session(subj_path)
            if not ses:
                print(f"  WARNING: No session found for {subj}")
                continue
            
            tract_dir = os.path.join(subj_path, ses, 'dwi', 'wholebrain', 'passthrough_tracts')
            
            subj_data = {}
            segments = ctrl_segments if group == 'CTRL' else graft_segments
            tgt_list = targets[group]
            
            for seg in segments:
                seg_data = {}
                for tgt in tgt_list:
                    # Construct tract filename
                    if group == 'CTRL':
                        tck_name = f"isocortex_{hemi}_{seg}_to_{tgt}_{hemi}_cleaned.tck"
                    else:
                        tck_name = f"seg9_graft_{hemi}_{seg}_to_{tgt}_{hemi}_cleaned.tck"
                    
                    tck_path = os.path.join(tract_dir, tck_name)
                    count = get_tck_count(tck_path)
                    
                    # Normalize target name
                    tgt_norm = tgt.lower().replace('_cortex', '')
                    if tgt_norm == 'palecortex':
                        tgt_norm = 'paleocortex'
                    seg_data[tgt_norm] = count
                subj_data[seg] = seg_data
            
            # Calculate total
            total = sum(sum(seg_data.values()) for seg_data in subj_data.values())
            subj_data['total'] = total
            
            data[group][subj.replace('sub-', '')] = subj_data
    
    return data

if __name__ == '__main__':
    print("Collecting RIGHT hemisphere data...")
    right_data = collect_hemisphere_data('right')
    print("\nRIGHT HEMISPHERE DATA:")
    print(json.dumps(right_data, indent=2))
    
    # Calculate totals
    print("\n=== RIGHT HEMISPHERE SUMMARY ===")
    for group in ['CTRL', 'AP', 'APTR']:
        totals = [v['total'] for v in right_data[group].values()]
        print(f"{group}: {totals} (sum={sum(totals)})")
