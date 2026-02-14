#!/usr/bin/env python3
"""
Example workflow: Draw ROIs in highres T2w and use for tractography

This demonstrates the complete workflow for using high-resolution anatomical
images to guide tractography in native DWI space.

Workflow:
1. Run two-stage coregistration (DWI → lowres → highres)
2. Draw ROIs in highres T2w space (best resolution)
3. Transform ROIs to DWI space using the transform chain
4. Use ROIs for tractography in native DWI space
"""

from pathlib import Path
from utils import DWIPreprocessor

# ============================================================================
# STEP 1: Run two-stage coregistration
# ============================================================================
print("=" * 80)
print("STEP 1: Two-stage coregistration")
print("=" * 80)
print("\nRun from command line:")
print("  python main.py --preprocess --steps coreg --first-subject\n")
print("This creates:")
print("  • DWI in lowres T2w space")
print("  • Transform chain: highres → lowres → DWI")
print("  • Transform file: *_highres2dwi.txt")
print()

# ============================================================================
# STEP 2: Draw ROIs in highres T2w space
# ============================================================================
print("=" * 80)
print("STEP 2: Draw ROIs in highres T2w (best resolution)")
print("=" * 80)
print("\nOption A - Using mrview (MRtrix3 viewer):")
print("  mrview data/bids/sub-XXX/ses-YYY/anat/*_acq-highres_T2w.nii.gz")
print("  • Tool → ROI editor")
print("  • Draw your ROI")
print("  • File → Save as... → roi_motor_cortex_highres.mif")
print()
print("Option B - Using ITK-SNAP or similar:")
print("  • Load: *_acq-highres_T2w.nii.gz")
print("  • Draw ROI and save as NIfTI")
print()

# ============================================================================
# STEP 3: Transform ROI from highres to DWI space
# ============================================================================
print("=" * 80)
print("STEP 3: Transform ROI to DWI space")
print("=" * 80)
print("\nExample code:")
print("""
from utils import DWIPreprocessor

# Initialize preprocessor
preprocessor = DWIPreprocessor()

# Transform ROI: highres → lowres → DWI
roi_dwi = preprocessor.transform_roi_to_dwi(
    subject_id='596-C-OO-NmM',
    session_id='240911103349',
    roi_highres=Path('roi_motor_cortex_highres.mif')
)

print(f"ROI transformed to DWI space: {roi_dwi}")
""")
print()

# ============================================================================
# STEP 4: Use ROI for tractography
# ============================================================================
print("=" * 80)
print("STEP 4: Use ROI for tractography")
print("=" * 80)
print("\nExample tractography command:")
print("""
tckgen \\
  data/derivatives/mrtrix_preproc/sub-XXX/ses-YYY/dwi/*_fod_wm.mif \\
  tracks_from_roi.tck \\
  -seed_image roi_motor_cortex_dwi.mif \\
  -mask *_mask.mif \\
  -select 10000 \\
  -angle 60 \\
  -step 0.02 \\
  -minlength 1 \\
  -maxlength 50
""")
print()

# ============================================================================
# COMPLETE EXAMPLE
# ============================================================================
print("=" * 80)
print("COMPLETE EXAMPLE (uncomment to run)")
print("=" * 80)
print("""
# Run this after:
# 1. Coregistration is complete
# 2. You have drawn an ROI in highres space

from pathlib import Path
from utils import DWIPreprocessor

# Configuration
SUBJECT_ID = '596-C-OO-NmM'
SESSION_ID = '240911103349'
ROI_HIGHRES = Path('my_roi_highres.mif')  # Your ROI drawn in highres

# Initialize
preprocessor = DWIPreprocessor()

# Transform ROI to DWI space
roi_dwi = preprocessor.transform_roi_to_dwi(
    subject_id=SUBJECT_ID,
    session_id=SESSION_ID,
    roi_highres=ROI_HIGHRES
)

print(f"✓ ROI ready for tractography: {roi_dwi}")

# Now use roi_dwi for tractography:
# python main.py --tracking --steps tckgen
# (modify tckgen step to use roi_dwi as seed/include/etc.)
""")
print()

# ============================================================================
# NOTES
# ============================================================================
print("=" * 80)
print("IMPORTANT NOTES")
print("=" * 80)
print("""
Registration Quality Check:
  • Overlay DWI on T2w: mrview *_t2w_lowres.mif -overlay.load *_coreg.mif
  • Check alignment at tissue boundaries
  • Verify no major shifts or rotations

ROI Drawing Tips:
  • Draw in highres for best precision (100µm vs 200µm)
  • Use anatomical landmarks visible in T2w
  • Multiple ROIs can be combined (seed + include/exclude)
  • ROI labels are preserved during transformation

Transform Chain:
  • highres → lowres: anatomical registration (same modality)
  • lowres → DWI: functional registration (b0 to T2w)
  • Chain is composed into single highres → DWI transform
  • Nearest-neighbor interpolation preserves ROI labels

Saved Files:
  • *_dwi2lowres.txt: DWI to lowres transform
  • *_lowres2highres.txt: lowres to highres transform
  • *_highres2dwi.txt: Combined transform chain (use this!)
  • *_t2w_lowres.mif: lowres T2w reference
  • *_t2w_highres.mif: highres T2w reference
  • *_coreg.mif: DWI in lowres space
""")
