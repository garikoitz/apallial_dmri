#!/usr/bin/env python3
"""Extract central AP slice from T2w with graft overlay, rotate 90° CW, save as sample.png."""
import subprocess
import numpy as np
import nibabel as nib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import tempfile, os

DERIV = "data/derivatives/mrtrix_preproc_flip1"
BASE = f"{DERIV}/sub-Apallial20R1F/ses-240909094232"
t2w = f"{BASE}/dwi/t2w_highres_flip.mif"
fod = f"{BASE}/dwi/wm_fod_norm_group.mif"
graft_l = f"{BASE}/manual_roi/fod_space/roi_Graft_left_flip_fod.mif"
graft_r = f"{BASE}/manual_roi/fod_space/roi_Graft_right_flip_fod.mif"
AP_SLICE = 14

def load_mif(path):
    """Convert .mif to nifti, load with nibabel, return numpy array."""
    tmp = tempfile.mktemp(suffix=".nii.gz")
    subprocess.run(["mrconvert", path, tmp, "-quiet", "-force"], check=True)
    img = nib.load(tmp)
    data = img.get_fdata()
    os.unlink(tmp)
    return data

# Regrid T2w to FOD space
t2w_fod_tmp = tempfile.mktemp(suffix=".mif")
subprocess.run(["mrgrid", t2w, "regrid", t2w_fod_tmp,
                "-template", fod, "-interp", "linear",
                "-quiet", "-force"], check=True)

# Load data
t2w_data = load_mif(t2w_fod_tmp)
os.unlink(t2w_fod_tmp)
graft_l_data = load_mif(graft_l)
graft_r_data = load_mif(graft_r)

print(f"T2w shape: {t2w_data.shape}")
print(f"Graft L shape: {graft_l_data.shape}")
print(f"AP slice: {AP_SLICE}")

# Extract AP slice (axis 0 = AP in flip1)
t2w_slice = t2w_data[AP_SLICE, :, :]    # DV × ML
graft_l_slice = graft_l_data[AP_SLICE, :, :]
graft_r_slice = graft_r_data[AP_SLICE, :, :]

# Combine grafts
graft_combined = np.zeros_like(graft_l_slice)
graft_combined[graft_l_slice > 0] = 1  # left
graft_combined[graft_r_slice > 0] = 2  # right

# Create figure
fig, ax = plt.subplots(1, 1, figsize=(8, 6))

# Show T2w (grayscale)
ax.imshow(t2w_slice.T, cmap='gray', origin='lower', aspect='auto',
          vmin=np.percentile(t2w_slice, 2), vmax=np.percentile(t2w_slice, 99))

# Overlay graft ROIs with transparency
graft_overlay = np.ma.masked_where(graft_combined.T == 0, graft_combined.T)
cmap_graft = ListedColormap(['#FF4444', '#4444FF'])  # red=left, blue=right
ax.imshow(graft_overlay, cmap=cmap_graft, origin='lower', aspect='auto',
          alpha=0.4, vmin=0.5, vmax=2.5)

ax.set_title(f'sub-Apallial20R1F — AP slice {AP_SLICE}\nGraft: red=left, blue=right',
             fontsize=12)
ax.axis('off')

plt.tight_layout()

# Save, then rotate 90° CW
tmp_path = "tmp/sample_raw.png"
fig.savefig(tmp_path, dpi=150, bbox_inches='tight', facecolor='black')
plt.close()

# Rotate 90° clockwise using PIL
from PIL import Image
img = Image.open(tmp_path)
img_rotated = img.rotate(-90, expand=True)
img_rotated.save("tmp/sample.png")
os.unlink(tmp_path)
print("Saved: tmp/sample.png")
