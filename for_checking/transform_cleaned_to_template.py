#!/usr/bin/env python3
"""
Transform cleaned tracts to group template space and open mrview for visualization.
"""

import subprocess
from pathlib import Path

# Base paths
DERIV = Path("/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc")
TEMPLATES = DERIV / "group_templates"

# Group definitions
GROUPS = {
    "CTRL": {
        "subjects": ["sub-129-1014-L1M", "sub-129-921-R1F", "sub-129-921-R1L1F", "sub-138-1010-L1F"],
        "tract_type": "Isocortex_Thalamus"
    },
    "AP": {
        "subjects": ["sub-596-C-OO-NmM", "sub-769-C-O5-R1M", "sub-769-C-OO-NmM", "sub-875-C-OO-NM"],
        "tract_type": "Graft_Thalamus"
    },
    "APTR": {
        "subjects": ["sub-AP72-862-901-NmF", "sub-AP72-862-901-R1M", "sub-Apallial17NmF", "sub-Apallial20R1F"],
        "tract_type": "Graft_Thalamus"
    }
}

# Session mapping
SESSIONS = {
    "sub-129-1014-L1M": "ses-260115113814",
    "sub-129-921-R1F": "ses-260115121231",
    "sub-129-921-R1L1F": "ses-260116085433",
    "sub-138-1010-L1F": "ses-260115101441",
    "sub-596-C-OO-NmM": "ses-240911103349",
    "sub-769-C-O5-R1M": "ses-240911115733",
    "sub-769-C-OO-NmM": "ses-240911111853",
    "sub-875-C-OO-NM": "ses-240911100253",
    "sub-AP72-862-901-NmF": "ses-251210093342",
    "sub-AP72-862-901-R1M": "ses-251210100043",
    "sub-Apallial17NmF": "ses-240816091633",
    "sub-Apallial20R1F": "ses-240909094232",
}

def transform_tract(tract_path, deform_path, output_path):
    """Transform tract to template space using tcktransform."""
    cmd = ["tcktransform", str(tract_path), str(deform_path), str(output_path), "-force"]
    print(f"  Transforming: {tract_path.name} -> {output_path.name}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    ERROR: {result.stderr}")
        return False
    return True

def main():
    for group_name, group_info in GROUPS.items():
        print(f"\n=== Processing {group_name} group ===")
        template_dir = TEMPLATES / group_name
        warped_dir = template_dir / "tracts_warped"
        warped_dir.mkdir(exist_ok=True)
        
        tract_type = group_info["tract_type"]
        
        for subject in group_info["subjects"]:
            session = SESSIONS.get(subject)
            if not session:
                print(f"  WARNING: No session found for {subject}")
                continue
            
            # Deformation field
            deform = warped_dir / f"{subject}_deformation.mif"
            if not deform.exists():
                print(f"  WARNING: No deformation field for {subject}")
                continue
            
            # Subject dwi directory
            dwi_dir = DERIV / subject / session / "dwi"
            
            for side in ["left", "right"]:
                # Source cleaned tract
                tract_name = f"tracks_{tract_type}_{side}_cleaned.tck"
                tract_path = dwi_dir / tract_name
                
                if not tract_path.exists():
                    print(f"  WARNING: {tract_name} not found for {subject}")
                    continue
                
                # Output warped tract
                output_name = f"{subject}_tracks_{tract_type}_{side}_cleaned_template.tck"
                output_path = warped_dir / output_name
                
                transform_tract(tract_path, deform, output_path)
        
        print(f"  Done with {group_name}")

if __name__ == "__main__":
    main()
