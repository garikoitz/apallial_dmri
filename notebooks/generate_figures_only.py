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
# tAP-hCO mouse, overlaid on a T2-weighted (T2W) image. Green streamlines
# indicate a dorso-ventral (DV) preference for water diffusion, red indicate an
# anterior-posterior (AP) preference, and blue indicate a medial-lateral (ML)
# preference. Dashed white lines delineate the location of the hCO graft.
#
# ### Figure 2h
# Top- Coronal T2W MRI scans in control, apallial, and tAP-hCO mice at 4-5
# months (organoid age post-differentiation and animal age) with fiber
# orientation distributions (FODs) estimated from superimposed diffusion-weighted
# imaging (DWI). FODs are coded with the same criteria as in panel g. The cortex
# or graft margin is indicated by a dashed line. The lack of an average diffusion
# direction in the apallial cavity is consistent with the absence of tissue and
# coherent microstructure. Scale bar, 1 mm. Bottom- Radar plots depict the
# principal directions of water diffusion for control, apallial, and tAP-hCO
# mice. Dashed lines delineate isocortex, apallial cavity, and graft in control,
# apallial, and tAP-hCO groups, respectively. Each plot depicts n= 4 mice per
# condition, differentiated by the shape of the points and shade of the
# respective color. See Extended Data Fig. 3h for the full legend.
#
# ### Extended Data Fig. 3i
# Coronal T2w MRI in tAP-hCO ~4.7 months post-differentiation with fiber
# orientation distributions (FODs) estimated from diffusion-weighted imaging
# (DWI) superimposed. Green FODs indicate a dorso-ventral (DV) preference for
# water diffusion, red: anterior-posterior (AP), and blue: medial-lateral (ML).
# Graft margin is indicated by white dashed line. Scale bar, 1mm.
#
# ### Extended Data Fig. 3j
# The distribution of dominant directions of water diffusion (see methods).
# rmANOVA: group main effect F (1, 6) = 3.000 (P=0.1340), direction main effect
# F (1.082, 6.494) = 6.347 (P=0.0409), interaction F (1.082, 6.494) = 7.399
# (P=0.0306); Holm-Sidak corrected posthoc tests for Control vs tAP-hCO: DV
# (P=0.0014), ML (P=0.5367), AP (P=0.1262).
#
# ### Extended Data Fig. 3k
# Radar plots depicting the principal directions of the fixels across 11
# consecutive slices. Rows represent the experimental groups: control (top),
# apallial (middle), and tAP-hCO (bottom). Columns represent the fixed, indexed
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
from pathlib import Path

# Resolve repository root by searching upward for config1.yaml.
_cwd = Path.cwd().resolve()
REPO_DIR = next((p for p in [_cwd, *_cwd.parents] if (p / "config1.yaml").exists()), None)
if REPO_DIR is None:
    raise RuntimeError("Could not find repository root containing config1.yaml")

sys.path.insert(0, str(REPO_DIR))
sys.path.insert(0, str(REPO_DIR / "functions"))

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
# tAP-hCO mouse, overlaid on a T2-weighted (T2W) image. Green streamlines
# indicate a dorso-ventral (DV) preference for water diffusion, red indicate an
# anterior-posterior (AP) preference, and blue indicate a medial-lateral (ML)
# preference. Dashed white lines delineate the location of the hCO graft.
#
# **Inputs**: Whole-brain tractogram (`.tck`), T2W background (`.mif`), graft
#   ROI for dashed outline
# **Output**: `reports/figures/figure_2g.png`

# %% [Figure 2g] Whole-brain tractogram with direction-colored streamlines
# TODO: implement Figure 2g here.
# Use FIGURES_DIR / "figure_2g.png" as the save path.
# Keep imports local to this cell so each figure is self-contained.
raise NotImplementedError("Figure 2g not yet implemented.")

# %% [markdown] Figure 2h — Coronal FODs + radar plots (top/bottom composite)
# ## Figure 2h — Coronal FODs + radar plots
#
# **Top**: Coronal T2W MRI scans in control, apallial, and tAP-hCO mice at 4-5
# months with fiber orientation distributions (FODs) estimated from DWI
# superimposed. FODs are color-coded: green = DV, red = AP, blue = ML. The
# cortex or graft margin is indicated by a dashed line. Scale bar, 1 mm.
#
# **Bottom**: Radar plots depicting the principal directions of water diffusion
# for control, apallial, and tAP-hCO mice. Dashed lines delineate isocortex,
# apallial cavity, and graft. Each plot depicts n=4 mice per condition,
# differentiated by point shape and color shade.
#
# **Inputs**: T2W images (`.mif`), FOD images (`.mif`), ROI outlines, per-subject
#   direction data (CSV or config-derived)
# **Output**: `reports/figures/figure_2h.png`

# %% [Figure 2h] Coronal FODs + radar plots composite
# TODO: implement Figure 2h here.
raise NotImplementedError("Figure 2h not yet implemented.")

# %% [markdown] Extended Data Fig. 3i — tAP-hCO coronal FODs with graft margin
# ## Extended Data Fig. 3i — tAP-hCO coronal FODs
#
# Coronal T2w MRI in tAP-hCO ~4.7 months post-differentiation with fiber
# orientation distributions (FODs) estimated from DWI superimposed. Green FODs
# indicate DV preference, red: AP, blue: ML. Graft margin is indicated by white
# dashed line. Scale bar, 1mm.
#
# **Inputs**: T2W image (`.mif`), FOD image (`.mif`), graft ROI outline
# **Output**: `reports/figures/figure_3i.png`

# %% [Extended Data Fig. 3i] tAP-hCO coronal FODs with graft margin
# TODO: implement Figure 3i here.
raise NotImplementedError("Figure 3i not yet implemented.")

# %% [markdown] Extended Data Fig. 3j — Dominant direction distribution (bar chart + stats)
# ## Extended Data Fig. 3j — Dominant direction distribution
#
# The distribution of dominant directions of water diffusion. rmANOVA: group
# main effect F(1,6)=3.000 (P=0.1340), direction main effect F(1.082,6.494)=6.347
# (P=0.0409), interaction F(1.082,6.494)=7.399 (P=0.0306). Holm-Sidak corrected
# posthoc tests for Control vs tAP-hCO: DV (P=0.0014), ML (P=0.5367), AP
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
# consecutive slices. Rows: control (top), apallial (middle), tAP-hCO (bottom).
# Columns: fixed, indexed slice position (Slices 1-11) across all subjects.
# Fixels were clustered into two main populations, aligned to a principal
# reference direction, weighted by peak amplitude, and normalized by total
# voxels. Each plot depicts n=4 mice per condition, differentiated by point
# shape and color shade.
#
# **Inputs**: Per-subject per-slice fixel direction data, group assignments
# **Output**: `reports/figures/figure_3k.png`

# %% [Extended Data Fig. 3k] Fixel radar plots across slices × groups
import sys, os, importlib
import numpy as np
from pathlib import Path

# ── Resolve derivatives ───────────────────────────────────────────────
# The config defaults to Linux paths; override with local newmask data.
DERIV_BASE = REPO_DIR / "data" / "derivatives" / "mrtrix_newmask"
if not DERIV_BASE.exists():
    raise RuntimeError(f"Derivatives not found: {DERIV_BASE}")

os.environ["MBRAIN_DERIVATIVES_DIR"] = str(DERIV_BASE)

import functions.fixel_roi_analysis as _fri
importlib.reload(_fri)
from functions.fixel_roi_analysis import (
    analyze_3group_kmeans,
    plot_3group_kmeans_comparison,
    plot_3group_component_lines,
)

print("Step 3k: 3-Group Fixel K-means Direction Comparison\n")

# ── Pass 1: compute data to find unified y_max ────────────────────────
combos = [("left", False), ("right", False)]
all_results = {}
global_max = 0.0

for hemi, use_shrink in combos:
    tag = f"{hemi}"
    r = analyze_3group_kmeans(
        config, hemisphere=hemi, use_shrink=False,
        output_dir=FIGURES_DIR, n_clusters=2, conda_env="mbrain",
        force=False, show=False, verbose=True, plot=False,
    )
    all_results[tag] = (hemi, False, r)
    if r:
        for group in ["CTRL", "AP", "APTR"]:
            for subj, data in r["all_group_data"].get(group, {}).items():
                for s_idx, res in data.items():
                    for sd in res.get("scaled_directions", []):
                        val = np.max(np.abs(sd))
                        if val > global_max:
                            global_max = val

y_max_unified = max(0.05, np.ceil(global_max * 20) / 20)
print(f"  Unified y_max = {y_max_unified:.4f} (global max = {global_max:.4f})\n")

# ── Pass 2: generate plots with forced y_max ──────────────────────────
for tag, (hemi, use_shrink, r) in all_results.items():
    print(f"\n{'=' * 60}")
    print(f"  {tag.upper()}")
    print(f"{'=' * 60}")
    if not r:
        print(f"  ⚠ No data for {tag}")
        continue

    tns = r["target_n_slices"]

    fig_path = FIGURES_DIR / f"figure_3k_radar_{hemi}.png"
    plot_3group_kmeans_comparison(
        r["all_group_data"], hemi, tns,
        use_shrink=False, save_path=fig_path, show=False,
        y_max=y_max_unified,
    )
    print(f"  ✅ Radar: {fig_path}")

    lines_path = FIGURES_DIR / f"figure_3k_components_{hemi}.png"
    plot_3group_component_lines(
        r["all_group_data"], hemi, tns,
        use_shrink=False, save_path=lines_path, show=False,
        y_max=y_max_unified,
    )
    print(f"  ✅ Lines: {lines_path}")

print("\nFigure 3k done.")

# %% [markdown] Summary
# ## Summary
#
# All figures are saved to `reports/figures/`. Re-run any individual cell to
# regenerate a single figure without touching the others.
    