# %% [markdown] Generate Paper Figures (Option A)
# # Generate Paper Figures (Option A)
#
# This notebook generates the figures used in the paper **directly from the
# already-installed data** (no preprocessing, no FOD, no tractography).
#
# It is the lightweight entry point for **Option A** in the README:
#
# > Install the pre-computed derivatives with `python utils/setup_data.py`,
# > then run only the figure-generation steps in
# > `notebooks/generate_figures_only.py`.
#
# Each cell below produces **one figure**. Run cells one at a time with
# `Shift+Enter` (in VS Code's notebook view) or run them all from top to
# bottom. Figures are saved into `reports/figures/` next to this notebook.
#
# ---
#
# **Preconditions**
#
# 1. The `mbrain` conda env is active (see `mbrain_env.yml`).
# 2. Data has been installed via `python utils/setup_data.py` (components
#    `local`, `pre_bids`, `pre_preproc`, `post`).
# 3. `config1.yaml` paths are correct for the current machine.
# 4. MRtrix3 3.0.4+ is available (required for any `.mif`/`.nii.gz` I/O).
#
# **Not required for this notebook**: ANTs, FreeSurfer, MATLAB, antspyx.
#
# ---
#
# ## Figures in the paper
#
# ### Figure 2g
# Whole-brain tractogram with 1M streamlines (10,000 shown) of a representative
# xenocortex (XCX) mouse, overlaid on a T2-weighted (T2W) image. Green streamlines
# indicate a dorso-ventral (DV) preference for water diffusion, red indicate an
# anterior-posterior (AP) preference, and blue indicate a medial-lateral (ML)
# preference. Dashed white lines delineate the location of the hCO graft.
#
# ### Figure 2h
# Top- Coronal T2W MRI scans in control, apallial, and xenocortex (XCX) mice at 4-5
# months (organoid age post-differentiation and animal age) with fiber
# orientation distributions (FODs) estimated from superimposed diffusion-weighted
# imaging (DWI). FODs are coded with the same criteria as in panel g. The cortex
# or graft margin is indicated by a dashed line. The lack of an average diffusion
# direction in the apallial cavity is consistent with the absence of tissue and
# coherent microstructure. Scale bar, 1 mm. Bottom- Radar plots depict the
# principal directions of water diffusion for control, apallial, and xenocortex (XCX)
# mice. Dashed lines delineate isocortex, apallial cavity, and graft in control,
# apallial, and xenocortex (XCX) groups, respectively. Each plot depicts n= 4 mice per
# condition, differentiated by the shape of the points and shade of the
# respective color. See Extended Data Fig. 3h for the full legend.
#
# ### Extended Data Fig. 3i
# Coronal T2w MRI in xenocortex (XCX) ~4.7 months post-differentiation with fiber
# orientation distributions (FODs) estimated from diffusion-weighted imaging
# (DWI) superimposed. Green FODs indicate a dorso-ventral (DV) preference for
# water diffusion, red: anterior-posterior (AP), and blue: medial-lateral (ML).
# Graft margin is indicated by white dashed line. Scale bar, 1mm.
#
# ### Extended Data Fig. 3j
# The distribution of dominant directions of water diffusion (see methods).
# rmANOVA: group main effect F (1, 6) = 3.000 (P=0.1340), direction main effect
# F (1.082, 6.494) = 6.347 (P=0.0409), interaction F (1.082, 6.494) = 7.399
# (P=0.0306); Holm-Sidak corrected posthoc tests for Control vs xenocortex (XCX): DV
# (P=0.0014), ML (P=0.5367), AP (P=0.1262).
#
# ### Extended Data Fig. 3k
# Radar plots depicting the principal directions of the fixels across 11
# consecutive slices. Rows represent the experimental groups: control (top),
# apallial (middle), and xenocortex (XCX) (bottom). Columns represent the fixed, indexed
# slice position (Slices 1-11) across all subjects. Note that the comparison in
# Fig. 2h displays a specific, per-subject anatomical landmark (the anterior
# margin of the thalamus), which maps variably to Slices 2-6 depending on the
# individual. This supplementary figure is not aligned to anatomical landmarks,
# but rather the relative position in the scan sequence. To quantify alignment,
# the fixels were clustered into two main populations, aligned to a principal
# reference direction, weighted by peak amplitude, and normalized by total
# voxels. Each plot depicts n=4 mice per condition, differentiated by the shape
# of the points, and shade of the respective color.
#
# ---

# %% [markdown] Setup — paths, config, output directory
# ## Setup
#
# Resolve the repository root, load `config1.yaml`, and prepare the output
# directory `reports/figures/`. This cell does not produce a figure.

# %% [Setup] Resolve repo, load config, prepare output dir
from __future__ import annotations

import sys
import os

from pathlib import Path


# Resolve repository root by searching upward for config1.yaml.
_cwd = Path.cwd().resolve()
REPO_DIR = next((p for p in [_cwd, *_cwd.parents] if (p / "config1.yaml").exists()), None)
if REPO_DIR is None:
    raise RuntimeError("Could not find repository root containing config1.yaml")

sys.path.insert(0, str(REPO_DIR))
sys.path.insert(0, str(REPO_DIR / "functions"))

# Override derivatives to local newmask data BEFORE loading config,
# so config['paths']['derivatives_flip1'] points to mrtrix_newmask.
_NEWMASK = REPO_DIR / "data" / "derivatives" / "mrtrix_newmask"
if not _NEWMASK.exists():
    raise RuntimeError(f"Derivatives not found: {_NEWMASK}")
os.environ["MBRAIN_DERIVATIVES_DIR"] = str(_NEWMASK)

from functions.config_loader import load_config

config = load_config(str(REPO_DIR / "config1.yaml"))

FIGURES_DIR = REPO_DIR / "reports" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

print(f"REPO_DIR    : {REPO_DIR}")
print(f"FIGURES_DIR : {FIGURES_DIR}")
print(f"Config keys : {list(config.keys())[:6]} ...")

# %% [markdown] Figure 2g — Whole-brain tractogram with direction-colored streamlines
# ## Figure 2g — Whole-brain tractogram (direction-colored)
#
# Whole-brain tractogram with 1M streamlines (10,000 shown) of a representative
# xenocortex (XCX) mouse, overlaid on a T2-weighted (T2W) image. Green streamlines
# indicate a dorso-ventral (DV) preference for water diffusion, red indicate an
# anterior-posterior (AP) preference, and blue indicate a medial-lateral (ML)
# preference. Dashed white lines delineate the location of the hCO graft.
#
# **Inputs**: Whole-brain tractogram (`.tck`), T2W background (`.mif`), graft
#   ROI for dashed outline
# **Output**: `reports/figures/figure_2g.png`

# %% [Figure 2g] Whole-brain tractogram with direction-colored streamlines
# Use FIGURES_DIR / "figure_2g.png" as the save path.
# Keep imports local to this cell so each figure is self-contained.

cmd = ["freeview",
  "/Users/glerma/soft/apallial_dmri/data/derivatives/mrtrix_newmask/sub-AP72-862-901-NmF/ses-251210093342/dwi/t2w_highres_flip.nii.gz",
  "-t", "/Users/glerma/soft/apallial_dmri/data/derivatives/mrtrix_newmask/sub-AP72-862-901-NmF/ses-251210093342/dwi/wholebrain/wholebrain_10k.trk"
]
print("Run the following command in a terminal to view the tractogram:")
print(" ".join(cmd))

# %% [markdown] Figure 2h — Coronal FODs + radar plots (top/bottom composite)
# ## Figure 2h — Coronal FODs + radar plots
#
# **Top**: Coronal T2W MRI scans in control, apallial, and xenocortex (XCX) mice at 4-5
# months with fiber orientation distributions (FODs) estimated from DWI
# superimposed. FODs are color-coded: green = DV, red = AP, blue = ML. The
# cortex or graft margin is indicated by a dashed line. Scale bar, 1 mm.
#
# **Bottom**: Radar plots depicting the principal directions of water diffusion
# for control, apallial, and xenocortex (XCX) mice. Dashed lines delineate isocortex,
# apallial cavity, and graft. Each plot depicts n=4 mice per condition,
# differentiated by point shape and color shade.
#
# **Inputs**: T2W images (`.mif`), FOD images (`.mif`), ROI outlines, per-subject
#   direction data (CSV or config-derived)
# **Output**: `reports/figures/figure_2h.png`

# %% [Figure 2h] Coronal FODs + radar plots composite
# Adapted from figs/make_4x3_figure_v09.py: 2-row × 3-column composite.
# Row 1 (top): FOD screenshots (coronal T2W + FOD glyphs, rotated 90° CW)
#             with dashed ROI margin overlay and 1 mm scale bar.
# Row 2 (bottom): Per-mouse radar plots (n=4 mice/group, bilateral L|R ROI,
#                 first-thalamus slice, K-means 2 clusters, DV-dominant solid).
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path

from functions.fixel_roi_analysis import (
    load_fixel_dir, load_roi_mask, extract_slicewise_vectors,
    slicewise_kmeans,
)
from notebooks.figure_helpers import (
    load_rotated, thin_radial_grid, draw_radar_mice,
    add_scale_bar, add_dashed_margin,
)

# --- Paths ---
TMP_DIR = REPO_DIR / "tmp"
DERIV = Path(config["paths"]["derivatives_flip1"])

# --- Representative subjects for FOD screenshots ---
REPR_SUBJECTS_2H = {
    "CTRL": {"id": "sub-138-1010-L1F", "ses": "ses-260115101441",
             "fod_slice": 14, "roi_name": "isocortex"},
    "AP":   {"id": "sub-769-C-O5-R1M",  "ses": "ses-240911115733",
             "fod_slice": 15, "roi_name": "graft"},
    "APTR": {"id": "sub-AP72-862-901-NmF", "ses": "ses-251210093342",
             "fod_slice": 14, "roi_name": "graft"},
}

# --- All subjects for radar (4 per group, bilateral ROI) ---
ALL_SUBJECTS_2H = {
    "CTRL": [
        {"id": "sub-129-921-R1F",   "ses": "ses-260115121231", "first_thal_fod": 10,
         "roi_l": "atlas/rois/fod_space/roi_isocortex_left.nii.gz",
         "roi_r": "atlas/rois/fod_space/roi_isocortex_right.nii.gz"},
        {"id": "sub-129-921-R1L1F", "ses": "ses-260116085433", "first_thal_fod": 10,
         "roi_l": "atlas/rois/fod_space/roi_isocortex_left.nii.gz",
         "roi_r": "atlas/rois/fod_space/roi_isocortex_right.nii.gz"},
        {"id": "sub-129-1014-L1M",  "ses": "ses-260115113814", "first_thal_fod": 9,
         "roi_l": "atlas/rois/fod_space/roi_isocortex_left.nii.gz",
         "roi_r": "atlas/rois/fod_space/roi_isocortex_right.nii.gz"},
        {"id": "sub-138-1010-L1F",  "ses": "ses-260115101441", "first_thal_fod": 9,
         "roi_l": "atlas/rois/fod_space/roi_isocortex_left.nii.gz",
         "roi_r": "atlas/rois/fod_space/roi_isocortex_right.nii.gz"},
    ],
    "AP": [
        {"id": "sub-596-C-OO-NmM", "ses": "ses-240911103349", "first_thal_fod": 13,
         "roi_l": "manual_roi/fod_space/roi_Graft_left_flip_fod.mif",
         "roi_r": "manual_roi/fod_space/roi_Graft_right_flip_fod.mif"},
        {"id": "sub-769-C-O5-R1M", "ses": "ses-240911115733", "first_thal_fod": 12,
         "roi_l": "manual_roi/fod_space/roi_Graft_left_flip_fod.mif",
         "roi_r": "manual_roi/fod_space/roi_Graft_right_flip_fod.mif"},
        {"id": "sub-769-C-OO-NmM", "ses": "ses-240911111853", "first_thal_fod": 11,
         "roi_l": "manual_roi/fod_space/roi_Graft_left_flip_fod.mif",
         "roi_r": "manual_roi/fod_space/roi_Graft_right_flip_fod.mif"},
        {"id": "sub-875-C-OO-NM",  "ses": "ses-240911100253", "first_thal_fod": 12,
         "roi_l": "manual_roi/fod_space/roi_Graft_left_flip_fod.mif",
         "roi_r": "manual_roi/fod_space/roi_Graft_right_flip_fod.mif"},
    ],
    "APTR": [
        {"id": "sub-AP72-862-901-NmF", "ses": "ses-251210093342", "first_thal_fod": 10,
         "roi_l": "manual_roi/fod_space/roi_Graft_left_flip_fod.mif",
         "roi_r": "manual_roi/fod_space/roi_Graft_right_flip_fod.mif"},
        {"id": "sub-AP72-862-901-R1M", "ses": "ses-251210100043", "first_thal_fod": 10,
         "roi_l": "manual_roi/fod_space/roi_Graft_left_flip_fod.mif",
         "roi_r": "manual_roi/fod_space/roi_Graft_right_flip_fod.mif"},
        {"id": "sub-Apallial17NmF",    "ses": "ses-240816091633", "first_thal_fod": 10,
         "roi_l": "manual_roi/fod_space/roi_Graft_left_flip_fod.mif",
         "roi_r": "manual_roi/fod_space/roi_Graft_right_flip_fod.mif"},
        {"id": "sub-Apallial20R1F",    "ses": "ses-240909094232", "first_thal_fod": 11,
         "roi_l": "manual_roi/fod_space/roi_Graft_left_flip_fod.mif",
         "roi_r": "manual_roi/fod_space/roi_Graft_right_flip_fod.mif"},
    ],
}

# --- Screenshot file map (user-captured PNGs) ---
# Run the mrview commands in Step 0 below to generate these files.
SCREENSHOT_MAP_2H = {
    ("CTRL", "fod"): REPO_DIR / "screenshots" / "ctrl_fod_v08_0000.png",
    ("AP", "fod"):   REPO_DIR / "screenshots" / "ap_fod_v08_0000.png",
    ("APTR", "fod"): REPO_DIR / "screenshots" / "aptr_fod_v08_0000.png",
}

# --- Per-group mouse color shades (4 shades each) ---
GROUP_MOUSE_COLORS_2H = {
    "CTRL": ["#1a1a1a", "#3d3d3d", "#5c5c5c", "#7a7a7a"],
    "AP":   ["#600030", "#800040", "#a0005c", "#c03080"],
    "APTR": ["#1a3fa0", "#3366cc", "#5588dd", "#88aaee"],
}
MOUSE_MARKERS_2H = ["o", "s", "^", "D"]

GROUPS_ORDER_2H = ["CTRL", "AP", "APTR"]
GROUP_TITLES_2H = ["Control (isocortex)", "Apallial (graft)", r"$t^{AP}$-hCO (graft + Tx)"]
GROUP_DISPLAY_2H = {"CTRL": "Control", "AP": "Apallial", "APTR": r"$t^{AP}$-hCO"}


# ============================================================
# Step 0: mrview screenshot capture instructions
# ============================================================
# Run these commands in a terminal to capture coronal FOD screenshots.
# Each opens mrview at the first-thalamus slice for the representative subject,
# with FOD glyphs overlaid on T2w. The -capture.grab flag auto-saves a PNG.
# After capturing, re-run Step 1–4 below to compose the figure.
#
# T2w voxel coordinates → FOD slice mapping:
#   CTRL: (24, 50, 53) → FOD slice 14  (sub-138-1010-L1F)
#   AP:   (25, 49, 60) → FOD slice 15  (sub-769-C-O5-R1M)
#   APTR: (24, 50, 54) → FOD slice 14  (sub-AP72-862-901-NmF)

SCREENSHOTS_DIR = REPO_DIR / "screenshots"
SCREENSHOTS_DIR.mkdir(exist_ok=True)

print("=" * 72)
print("STEP 0: Capture FOD screenshots with mrview")
print("=" * 72)
print(f"Screenshots will be saved to: {SCREENSHOTS_DIR}")
print()
print("Run these commands in a terminal (micromamba activate mbrain first):")
print()

for group in GROUPS_ORDER_2H:
    info = REPR_SUBJECTS_2H[group]
    sid = info["id"]
    ses = info["ses"]
    base = DERIV / sid / ses / "dwi"
    voxel_str = {
        "CTRL": "24,50,53",
        "AP":   "25,49,60",
        "APTR": "24,50,54",
    }[group]
    prefix = {"CTRL": "ctrl", "AP": "ap", "APTR": "aptr"}[group]

    cmd = (
        f'mrview "{base}/t2w_highres_flip.mif" \\\n'
        f'  -odf.load_sh "{base}/wm_fod_norm_group.mif" \\\n'
        f'  -voxel {voxel_str} \\\n'
        f'  -plane 0 \\\n'
        f'  -noannotations \\\n'
        f'  -capture.folder "{SCREENSHOTS_DIR}" \\\n'
        f'  -capture.prefix {prefix}_fod_v08_ \\\n'
        f'  -capture.grab \\\n'
        f'  -exit'
    )
    print(f"  # {group} — {sid} (FOD slice {info['fod_slice']})")
    print(f"  {cmd}")
    print()

print("After capturing, re-run the cell below to compose the figure.")
print()

# %% [markdown] Figure 2h (bottom row) — Compose 2×3 figure from screenshots + radar
# ## Figure 2h — Compose figure (bottom row + SVG/PNG)
#
# Run this cell **after** capturing the mrview screenshots in the cell above.
# It loads the screenshot PNGs, generates per-mouse radar data from fixels,
# composes the 2-row × 3-column figure, and saves both PNG and SVG.

# %% [Figure 2h — compose] Load screenshots + radar data + compose + save
# Steps 1–4: load FOD screenshots, generate per-mouse radar data (bilateral ROI),
#            compute data-driven y_max, compose 2×3 figure, save PNG + SVG.

# ============================================================
# Step 1: Load FOD screenshots
# ============================================================
print("Loading FOD screenshots...")
screenshots_2h = {}
for key, path in SCREENSHOT_MAP_2H.items():
    if path.exists():
        screenshots_2h[key] = load_rotated(path)
        print(f"  ✓ {key}: {path.name}  shape={screenshots_2h[key].shape}")
    else:
        print(f"  ✗ {key}: {path.name} NOT FOUND")
        screenshots_2h[key] = np.zeros((400, 400, 3), dtype=np.uint8)

# ============================================================
# Step 2: Generate per-mouse radar data (bilateral ROI)
# ============================================================
print("\n--- Generating per-mouse radar data (bilateral ROI) ---")
per_mouse_data_2h = {}

for group in GROUPS_ORDER_2H:
    mice_results = []
    print(f"\n  {group}:")
    for info in ALL_SUBJECTS_2H[group]:
        sid = info["id"]
        ses = info["ses"]
        fod_slice = info["first_thal_fod"]
        subj_dir = DERIV / sid / ses
        fixel_dir = subj_dir / "dwi" / "fixels_group"
        roi_l = subj_dir / info["roi_l"]
        roi_r = subj_dir / info["roi_r"]

        if not (fixel_dir / "index.mif").exists():
            print(f"    ⚠ Skipping {sid}: no fixels")
            continue
        if not roi_l.exists() or not roi_r.exists():
            print(f"    ⚠ Skipping {sid}: ROI not found")
            continue

        print(f"    {sid}: slice={fod_slice}...")
        fixel_data = load_fixel_dir(fixel_dir)
        mask_l = load_roi_mask(roi_l)
        mask_r = load_roi_mask(roi_r)
        roi_mask = mask_l | mask_r
        n_total = int(roi_mask.sum())
        print(f"      ROI voxels L|R={n_total}")

        slice_vecs, slice_peaks, slice_voxels = extract_slicewise_vectors(
            fixel_data, roi_mask, slice_axis=0)

        if fod_slice not in slice_vecs:
            avail = sorted(slice_vecs.keys())
            if not avail:
                print(f"      ⚠ No slices with fixels!")
                continue
            nearest = min(avail, key=lambda s: abs(s - fod_slice))
            print(f"      ⚠ Slice {fod_slice} → {nearest}")
            fod_slice = nearest

        vecs = slice_vecs[fod_slice]
        peaks = slice_peaks[fod_slice]
        n_vox = slice_voxels.get(fod_slice, len(vecs))

        if len(vecs) < 4:
            print(f"      ⚠ Only {len(vecs)} fixels, skipping")
            continue

        km_results = slicewise_kmeans(
            {fod_slice: vecs}, n_clusters=2,
            slice_afds={fod_slice: peaks},
            slice_roi_voxels={fod_slice: n_vox},
        )
        km_slice = km_results[fod_slice]

        dv_vals = [np.abs(sd[1]) for sd in km_slice["scaled_directions"]]
        best_k = int(np.argmax(dv_vals))
        # other_k = 1 - best_k
        scaled_dir = km_slice["scaled_directions"][best_k]
        # scaled_dir_other = km_slice["scaled_directions"][other_k]
        n_cluster = km_slice["cluster_sizes"][best_k]

        short = sid.replace("sub-", "")
        mice_results.append({
            "label": short,
            "scaled_dir": scaled_dir,
            # "scaled_dir_other": scaled_dir_other,
            "n_fixels": n_cluster,
        })
        print(f"      DV-dom={scaled_dir}, n={n_cluster}")

    per_mouse_data_2h[group] = mice_results

# ============================================================
# Step 3: Compute data-driven y_max
# ============================================================
all_maxes_2h = []
for group, mice in per_mouse_data_2h.items():
    for md in mice:
        all_maxes_2h.append(np.max(np.abs(md["scaled_dir"])))
        # if md.get("scaled_dir_other") is not None:
        #     all_maxes_2h.append(np.max(np.abs(md["scaled_dir_other"])))
y_max_2h = max(all_maxes_2h) if all_maxes_2h else 1.0
print(f"\nRadar y_max: {y_max_2h:.4f}")

# ============================================================
# Step 4: Compose 2×3 figure
# ============================================================
print("\n--- Composing 2×3 figure ---")
fig_2h = plt.figure(figsize=(16, 12))
gs_2h = fig_2h.add_gridspec(2, 3, hspace=0.25, wspace=0.15,
                            height_ratios=[1.0, 0.95],
                            top=0.92, bottom=0.04, left=0.06, right=0.96)

# Row 1: FOD screenshots with dashed margin + scale bar
for col, group in enumerate(GROUPS_ORDER_2H):
    ax = fig_2h.add_subplot(gs_2h[0, col])
    img = screenshots_2h.get((group, "fod"), np.zeros((10, 10, 3), dtype=np.uint8))
    ax.imshow(img)
    repr_info = REPR_SUBJECTS_2H[group]
    ax.set_title(f"{GROUP_TITLES_2H[col]}\n{repr_info['id'].replace('sub-', '')}",
                 fontsize=11, fontweight="bold", pad=8)
    ax.axis("off")
    # Dashed margin overlay (white dashed rectangle inset)
    add_dashed_margin(ax, color="white", linewidth=1.5, linestyle="--", inset_frac=0.06)
    # Scale bar (1 mm, bottom-right, white)
    add_scale_bar(ax, bar_mm=1.0, px_per_mm=47.0, position="br",
                  color="white", linewidth=2.5, pad_frac=0.05)
    if col == 0:
        ax.set_ylabel("FODs", fontsize=13, fontweight="bold", labelpad=15)

# Row 2: Per-mouse radar (4 mice, bilateral ROI)
for col, group in enumerate(GROUPS_ORDER_2H):
    ax = fig_2h.add_subplot(gs_2h[1, col], projection="polar")
    mice = per_mouse_data_2h.get(group, [])
    title = GROUP_DISPLAY_2H[group]
    if mice:
        draw_radar_mice(ax, mice, title=title, y_max=y_max_2h, group=group,
                        group_mouse_colors=GROUP_MOUSE_COLORS_2H,
                        mouse_markers=MOUSE_MARKERS_2H)
    else:
        ax.set_title(f"{group}: no data", fontsize=10)

# Row labels
fig_2h.text(0.015, 0.72, "FOD Glyphs\n(coronal)", fontsize=11, fontweight="bold",
            rotation=90, va="center")
fig_2h.text(0.015, 0.30, "Individual Mice\n(first-thalamus,\nbilateral)",
            fontsize=10, fontweight="bold", rotation=90, va="center")

fig_2h.suptitle("FOD Direction Analysis — Control vs Apallial vs "
                r"$t^{AP}$-hCO (bilateral)", fontsize=14, fontweight="bold", y=0.97)

# Save PNG
out_png_2h = FIGURES_DIR / "figure_2h.png"
fig_2h.savefig(out_png_2h, dpi=150, bbox_inches="tight", facecolor="white")
print(f"\n✅ Saved: {out_png_2h}")

# Save SVG (try cairo, fallback to fonttype=none)
out_svg_2h = FIGURES_DIR / "figure_2h.svg"
try:
    import cairo
    from matplotlib.backends.backend_cairo import FigureCanvasCairo
    canvas_orig = fig_2h.canvas
    FigureCanvasCairo(fig_2h)
    fig_2h.savefig(out_svg_2h, format="svg", bbox_inches="tight", facecolor="white")
    fig_2h.set_canvas(canvas_orig)
    print(f"✅ Saved (cairo SVG): {out_svg_2h}")
except Exception as e:
    print(f"⚠️  Cairo SVG failed ({e}), using svg.fonttype=none")
    with matplotlib.rc_context({"svg.fonttype": "none"}):
        fig_2h.savefig(out_svg_2h, format="svg", bbox_inches="tight", facecolor="white")
    print(f"✅ Saved (SVG fonttype=none): {out_svg_2h}")

plt.close(fig_2h)
print("Figure 2h complete.")

# %% [markdown] Extended Data Fig. 3i — xenocortex (XCX) coronal FODs with graft margin
# ## Extended Data Fig. 3i — xenocortex (XCX) coronal FODs
#
# Coronal T2w MRI in xenocortex (XCX) ~4.7 months post-differentiation with fiber
# orientation distributions (FODs) estimated from DWI superimposed. Green FODs
# indicate DV preference, red: AP, blue: ML. Graft margin is indicated by white
# dashed line. Scale bar, 1mm.
#
# **Inputs**: T2W image (`.mif`), FOD image (`.mif`), graft ROI outline
# **Output**: `reports/figures/figure_3i.png`

# %% [Extended Data Fig. 3i] xenocortex (XCX) coronal FODs with graft margin
# TODO: implement Figure 3i here.
raise NotImplementedError("Figure 3i not yet implemented.")

# %% [markdown] Extended Data Fig. 3j — Dominant direction distribution (bar chart + stats)
# ## Extended Data Fig. 3j — Dominant direction distribution
#
# The distribution of dominant directions of water diffusion. rmANOVA: group
# main effect F(1,6)=3.000 (P=0.1340), direction main effect F(1.082,6.494)=6.347
# (P=0.0409), interaction F(1.082,6.494)=7.399 (P=0.0306). Holm-Sidak corrected
# posthoc tests for Control vs xenocortex (XCX): DV (P=0.0014), ML (P=0.5367), AP
# (P=0.1262).
#
# **Inputs**: Per-subject dominant direction data, group assignments from config
# **Output**: `reports/figures/figure_3j.png`

# %% [Extended Data Fig. 3j] Dominant direction distribution with stats
# TODO: implement Figure 3j here.
raise NotImplementedError("Figure 3j not yet implemented.")

# %% [markdown] Extended Data Fig. 3k — Fixel radar plots across 11 slices × 3 groups
# ## Extended Data Fig. 3k — Fixel radar plots (11 slices × 3 groups)
#
# Radar plots depicting the principal directions of the fixels across 11
# consecutive slices. Rows: control (top), apallial (middle), xenocortex (XCX) (bottom).
# Columns: fixed, indexed slice position (Slices 1-11) across all subjects.
# Fixels were clustered into two main populations, aligned to a principal
# reference direction, weighted by peak amplitude, and normalized by total
# voxels. Each plot depicts n=4 mice per condition, differentiated by point
# shape and color shade.
#
# **Inputs**: Per-subject per-slice fixel direction data, group assignments
# **Output**: `reports/figures/figure_3k.png`

# %% [Extended Data Fig. 3k] Fixel radar plots across slices × groups
# Adapted from figs/make_3group_v11.py: bilateral (L|R) ROIs, newmask pipeline,
# per-group color shades (grey=CTRL, magenta=AP, blue=APTR), per-mouse markers,
# white radar backgrounds, thin radial grid, slice 0 dropped, Dir1 only.
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path

from functions.fixel_roi_analysis import (
    load_fixel_dir, load_roi_mask, extract_slicewise_vectors,
    slicewise_kmeans,
)
from notebooks.figure_helpers import (
    get_roi_paths, missing_fixel_inputs, validate_figure_3k_inputs,
    thin_radial_grid, load_all_group_data, compute_global_max,
    plot_3group_radar_v11,
)

print("Step 3k: 3-Group Fixel K-means Direction Comparison (v11, bilateral)\n")

# ── Derivatives root (newmask pipeline) ───────────────────────────────
DERIV = Path(config['paths']['derivatives_flip1'])
if not DERIV.exists():
    raise RuntimeError(f"Derivatives not found: {DERIV}")

# ── Subject definitions (from config1.yaml groups) ───────────────────
SUBJECT_MAP = {
    'CTRL': [
        ('sub-129-921-R1F',   'ses-260115121231'),
        ('sub-129-921-R1L1F', 'ses-260116085433'),
        ('sub-129-1014-L1M',  'ses-260115113814'),
        ('sub-138-1010-L1F',  'ses-260115101441'),
    ],
    'AP': [
        ('sub-596-C-OO-NmM', 'ses-240911103349'),
        ('sub-769-C-O5-R1M', 'ses-240911115733'),
        ('sub-769-C-OO-NmM', 'ses-240911111853'),
        ('sub-875-C-OO-NM',  'ses-240911100253'),
    ],
    'APTR': [
        ('sub-AP72-862-901-NmF', 'ses-251210093342'),
        ('sub-AP72-862-901-R1M', 'ses-251210100043'),
        ('sub-Apallial17NmF',    'ses-240816091633'),
        ('sub-Apallial20R1F',    'ses-240909094232'),
    ],
}

# Per-group mouse color shades (matching 4x3 v09 style)
GROUP_MOUSE_COLORS = {
    'CTRL': ['#1a1a1a', '#3d3d3d', '#5c5c5c', '#7a7a7a'],   # dark grey shades
    'AP':   ['#600030', '#800040', '#a0005c', '#c03080'],     # magenta shades
    'APTR': ['#1a3fa0', '#3366cc', '#5588dd', '#88aaee'],     # blue shades
}
MOUSE_MARKERS = ['o', 's', '^', 'D']


# ============================================================
# Helpers — imported from notebooks/figure_helpers.py
# (get_roi_paths, missing_fixel_inputs, validate_figure_3k_inputs,
#  thin_radial_grid, load_all_group_data, compute_global_max,
#  plot_3group_radar_v11)
# ============================================================


# ── Load all data (bilateral ROIs, slice 0 dropped) ──────────────────
validate_figure_3k_inputs(DERIV, SUBJECT_MAP)
all_group_data, target_n_slices = load_all_group_data(DERIV, SUBJECT_MAP, verbose=True)
y_max = compute_global_max(all_group_data)

plot_3group_radar_v11(
    all_group_data, target_n_slices, y_max,
    group_mouse_colors=GROUP_MOUSE_COLORS,
    mouse_markers=MOUSE_MARKERS,
    save_png=FIG_DIR / 'fig3k_3group_radar_v11.png',
    save_svg=FIG_DIR / 'fig3k_3group_radar_v11.svg',
)
