# %% [markdown] MBRAIN Pipeline: New Brain Mask for All Groups
# # MBRAIN Pipeline: New Brain Mask (Atlas-based) + Re-run Pipeline
#
# This notebook generates **new brain masks** using single-atlas segmentation
# (ANTs SyN registration) for ALL subjects across ALL groups (CTRL, AP, APTR),
# then re-runs the downstream pipeline steps that depend on the mask.
#
# **Workflow:**
# 1. Setup — load config, enumerate all subjects
# 2. Atlas-based mask — register atlas T2w+mask to each target
# 3. Visual QC — scroll through each mask overlaid on T2w
# 4. Manual Correction — fix masks for grafted mice (AP/APTR) where needed
# 5. Re-run downstream — bias correction, registration, FOD, tractography, etc.
#
# **Important:** Run each cell on lmc02. Requires modules: mrtrix3, ants, freesurfer, matlab.
#
# ---

# %% [markdown] Step 0: Setup & Configuration
# ## Step 0: Setup & Configuration
#
# Load config, enumerate all subjects from all groups.

# %% [Step 0] Setup — Load Config & Build Subject List
import os, sys, subprocess, glob, json
from pathlib import Path
from datetime import datetime

# ============================================================
# MACHINE-SPECIFIC: Set REPO_DIR for lmc02
# ============================================================
REPO_DIR = Path("/export/home/glerma/glerma/soft/mbrain")
DATA_BASE = Path("/export/home/glerma/glerma/PROYECTOS/mbrain")
MODULE_LOAD_FS = "module load freesurfer/7.4.1-Rocky8"
MODULE_LOAD_MRTRIX = "module load mrtrix3/3.0.3"
MODULE_LOAD_MATLAB = "module load matlab/R2024A"
MODULE_LOAD_ANTS = "module load ants/latest"
MODULES = f"{MODULE_LOAD_FS} && {MODULE_LOAD_MRTRIX} && {MODULE_LOAD_ANTS}"

sys.path.insert(0, str(REPO_DIR))
sys.path.insert(0, str(REPO_DIR / "functions"))

# Load config (now resolves {variable} placeholders automatically)
from functions.config_loader import load_config, get_subject_groups
cfg = load_config()
paths = cfg['paths']
groups = get_subject_groups()

# Config sub-sections used in downstream steps
tck = cfg['tractography']
wb  = cfg.get('wholebrain_tractography', {})
pt  = cfg.get('passthrough_extraction', {})
mat = cfg['matlab_cleaning']
sk  = cfg.get('shrink_analysis', {})

# For AP/APTR: segmentation mode (9 = 3 A-P x 3 D-V — recommended)
SEG_MODE = 9

# Source (original) and destination (new mask) derivative directories
DERIV_SRC = Path(paths['derivatives_flip1'])  # mrtrix_preproc_flip1
DERIV_BASE = DATA_BASE / "data" / "derivatives" / "mrtrix_newmask"

print(f"Source derivatives: {DERIV_SRC}")
print(f"New mask derivatives: {DERIV_BASE}")

# %% [Step 0b] Copy pre-mask files to mrtrix_newmask
# ============================================================
# Copy the full mrtrix_preproc_flip1 tree to mrtrix_newmask.
# This preserves the original pipeline outputs while we work
# on new masks in a separate folder.
#
# Uses rsync to be resumable and efficient (skips existing files).
# ============================================================

print(f"{'=' * 70}")
print(f" COPYING mrtrix_preproc_flip1 → mrtrix_newmask")
print(f"{'=' * 70}")
print(f"  Source: {DERIV_SRC}")
print(f"  Dest:   {DERIV_BASE}")

if DERIV_BASE.exists() and any(DERIV_BASE.iterdir()):
    print(f"\n  ⚠️  Destination already exists and is not empty.")
    print(f"      Rsync will only copy new/updated files (safe to re-run).")

os.makedirs(DERIV_BASE, exist_ok=True)

result = subprocess.run(
    f"rsync -av --progress {DERIV_SRC}/ {DERIV_BASE}/",
    shell=True, capture_output=True, text=True, timeout=3600
)

if result.returncode == 0:
    # Count subjects copied
    copied_subs = [d for d in DERIV_BASE.iterdir() if d.is_dir() and d.name.startswith("sub-")]
    print(f"\n✅ Copy complete: {len(copied_subs)} subject folders in {DERIV_BASE.name}")
else:
    print(f"\n❌ rsync failed (exit code {result.returncode})")
    print(result.stderr[:500])

# %% [Step 0c] Build subject list from mrtrix_newmask
# Build flat list of all subjects with their group, session, and paths
# All paths now point to the mrtrix_newmask directory
# Exclude test mice — only real subjects (4 per group = 12 total)
SKIP_IDS = {'sub-APtestmouseNmM'}

subjects = []
for group_name in ["CTRL", "AP", "APTR"]:
    for s in groups.get(group_name, []):
        sid = s['id']
        if sid in SKIP_IDS:
            continue
        ses = s['session']
        base_dir = DERIV_BASE / sid / ses
        dwi_dir = base_dir / "dwi"
        wb_dir = dwi_dir / "wholebrain"
        subjects.append({
            'id': sid,
            'session': ses,
            'group': group_name,
            'base': base_dir,
            'dwi_dir': dwi_dir,
            'wb_dir': wb_dir,
            't2w': dwi_dir / "t2w_highres_flip.mif",
            'mask_old': dwi_dir / "mask.mif",
            'mask_anat_old': dwi_dir / "mask_anat.mif",
            'mask_new': dwi_dir / "mask_new.nii.gz",
            'mask_new_mif': dwi_dir / "mask_new.mif",
            'mask_new_edited': dwi_dir / "mask_new_edited.mif",
            'wb_tck': wb_dir / "wholebrain_1000k.tck",
            'sift2_weights': wb_dir / "sift2_weights.txt",
            'sift2_mu': wb_dir / "sift2_mu.txt",
            'pt_dir': wb_dir / "passthrough_tracts",
            'pt_shrink_dir': wb_dir / "passthrough_tracts_shrink",
            'fixel_dir': dwi_dir / "fixels_group",
            'afd_voxel': dwi_dir / "afd_voxel_sum.mif",
        })

print(f"{'=' * 70}")
print(f" SUBJECTS TO PROCESS: {len(subjects)}")
print(f"{'=' * 70}")
for s in subjects:
    exists = "✅" if s['t2w'].exists() else "❌"
    print(f"  {exists} [{s['group']:>4}] {s['id']:<30} {s['session']}")
print(f"{'=' * 70}")
print(f"DERIV_BASE: {DERIV_BASE}")

# %% [Step 0] Helper Functions
def run(cmd, check=True, timeout=None, verbose=True):
    """Run a shell command, print it, and return the result."""
    if verbose:
        print(f"\n>>> {cmd}\n")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
    if result.stdout.strip() and verbose:
        print(result.stdout.strip())
    if result.stderr.strip() and verbose:
        print("STDERR:", result.stderr.strip())
    if check and result.returncode != 0:
        print(f"⚠️ Command failed with exit code {result.returncode}")
    return result

def check_file(path, label=""):
    p = Path(path)
    if p.exists():
        size = p.stat().st_size / (1024*1024)
        print(f"  ✅ {label or p.name}: {size:.1f} MB")
        return True
    else:
        print(f"  ❌ MISSING: {label or p.name}")
        return False

def launch_gui(cmd):
    """Launch a GUI command (mrview, etc.) in the background without blocking."""
    print(f">>> {cmd}")
    subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

# ── Parallel runner for lmc02 (40 cores) ──
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

MAX_WORKERS = min(12, multiprocessing.cpu_count())  # 12 subjects = 12 workers max

# ── Log directory for parallel runs ──
LOG_DIR = DERIV_BASE / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)

def _run_subject_cmd(args):
    """Worker function for parallel execution.
    Saves full stdout/stderr to LOG_DIR/<label>_<sid>.log
    Returns (subject_id, returncode, stdout_tail, stderr_tail)."""
    sid, cmd, timeout, log_dir, label = args
    try:
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
        # Write full log to file
        from datetime import datetime
        log_file = Path(log_dir) / f"{label}_{sid}.log"
        with open(log_file, 'w') as f:
            f.write(f"{'=' * 70}\n")
            f.write(f"Subject: {sid}\n")
            f.write(f"Label:   {label}\n")
            f.write(f"Time:    {datetime.now().isoformat()}\n")
            f.write(f"RC:      {r.returncode}\n")
            f.write(f"Command: {cmd[:500]}\n")
            f.write(f"{'=' * 70}\n\n")
            f.write("=== STDOUT ===\n")
            f.write(r.stdout)
            f.write("\n\n=== STDERR ===\n")
            f.write(r.stderr)
        return (sid, r.returncode,
                r.stdout[-500:] if len(r.stdout) > 500 else r.stdout,
                r.stderr[-500:] if len(r.stderr) > 500 else r.stderr)
    except subprocess.TimeoutExpired:
        return (sid, -1, '', 'TIMEOUT')
    except Exception as e:
        return (sid, -1, '', str(e))

def run_parallel(tasks, label="Processing"):
    """Run a list of (subject_id, cmd, timeout) in parallel.
    Logs full output per subject to LOG_DIR/<label>_<sid>.log.
    Returns dict with 'success' and 'failed' lists of subject IDs."""
    results = {'success': [], 'failed': []}
    # Augment tasks with log_dir and label for the worker
    augmented = [(sid, cmd, to, str(LOG_DIR), label.replace(' ', '_'))
                 for sid, cmd, to in tasks]
    print(f"  Launching {len(tasks)} tasks with {MAX_WORKERS} parallel workers...")
    print(f"  Logs → {LOG_DIR}/")
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(_run_subject_cmd, t): t[0] for t in augmented}
        for future in as_completed(futures):
            sid = futures[future]
            try:
                sid_r, rc, stdout, stderr = future.result()
                if rc == 0:
                    results['success'].append(sid_r)
                    print(f"  ✅ {sid_r}: done")
                else:
                    results['failed'].append(sid_r)
                    print(f"  ❌ {sid_r}: FAILED (rc={rc})")
                    if stderr.strip():
                        print(f"     {stderr.strip()[:200]}")
                    print(f"     Full log: {LOG_DIR}/{label.replace(' ', '_')}_{sid_r}.log")
            except Exception as e:
                results['failed'].append(sid)
                print(f"  ❌ {sid}: EXCEPTION {e}")
    print(f"\n  {label}: ✅ {len(results['success'])} / ❌ {len(results['failed'])}")
    return results

print("Helpers loaded ✅")

# %% [Step 0] Verify tools: MRtrix3, ANTs, FreeSurfer, MATLAB
result_fs = run(f"{MODULE_LOAD_FS} && mri_info --version", check=False, verbose=False)
if result_fs.returncode == 0:
    print(f"✅ FreeSurfer found (module: {MODULE_LOAD_FS})")
else:
    print(f"❌ FreeSurfer NOT found")
    print(f"   Module load: {MODULE_LOAD_FS}")
    print(f"   stderr: {result_fs.stderr[:200]}")

result_mr = run(f"{MODULE_LOAD_MRTRIX} && mrconvert --version", check=False, verbose=False)
if result_mr.returncode == 0:
    print(f"✅ MRtrix3 found (module: {MODULE_LOAD_MRTRIX})")
else:
    print(f"❌ MRtrix3 NOT found")
    print(f"   Module load: {MODULE_LOAD_MRTRIX}")
    print(f"   stderr: {result_mr.stderr[:200]}")

result_ants = run(f"{MODULE_LOAD_ANTS} && N4BiasFieldCorrection --version", check=False, verbose=False)
if result_ants.returncode == 0 or 'ANTs' in (result_ants.stdout + result_ants.stderr):
    print(f"✅ ANTs found (module: {MODULE_LOAD_ANTS})")
else:
    print(f"❌ ANTs NOT found")
    print(f"   Module load: {MODULE_LOAD_ANTS}")
    print(f"   stderr: {result_ants.stderr[:200]}")

result_mat = run(f"{MODULE_LOAD_MATLAB} && matlab -batch \"disp('MATLAB OK')\"", check=False, verbose=False)
if result_mat.returncode == 0 and 'MATLAB OK' in result_mat.stdout:
    print(f"✅ MATLAB found (module: {MODULE_LOAD_MATLAB})")
else:
    print(f"❌ MATLAB NOT found")
    print(f"   Module load: {MODULE_LOAD_MATLAB}")
    print(f"   stderr: {result_mat.stderr[:200]}")


# %% [markdown] Step 0.9: Delete existing mask_new files (re-run safety)
# ---
# ## Step 0.9: Delete existing mask_new files
#
# If re-running Step 1 with new parameters, first delete existing masks
# so they aren't skipped. Only deletes mask_new.* — NOT mask_new_edited.mif.

# %% [Step 0.9] Delete existing mask_new (but NOT mask_new_edited) to allow re-run
deleted = 0
for s in subjects:
    for f in [s['mask_new'], s['mask_new_mif']]:
        if f.exists():
            f.unlink()
            deleted += 1
            print(f"  🗑️ Deleted: {f.name} ({s['id']})")
print(f"\nDeleted {deleted} mask files. mask_new_edited.mif files preserved.")

# %% [markdown] Step 1: Generate Brain Masks (Group-specific atlas registration)
# ---
# ## Step 1: Generate Brain Masks — Group-Specific Atlas Segmentation
#
# Each group has its own **manually corrected atlas** subject whose mask is
# registered to the remaining subjects in that group via ANTs SyN:
#
# | Group | Atlas subject           | # targets |
# |-------|-------------------------|-----------|
# | CTRL  | sub-129-921-R1F         | 3         |
# | AP    | sub-596-C-OO-NmM       | 3         |
# | APTR  | sub-AP72-862-901-NmF   | 3         |
#
# Pipeline per target:
#   1. **N4BiasFieldCorrection** (ANTs) — corrects intensity inhomogeneity
#   2. **antsRegistrationSyN.sh** (ANTs) — rigid + affine + SyN nonlinear
#   3. **antsApplyTransforms** (ANTs) — warp atlas mask (NearestNeighbor)
#   4. **maskfilter** (MRtrix) — keep largest connected component
#
# Atlas subjects keep their own `mask_new_edited.mif` as-is.
#
# **Output:** `mask_new.mif` and `mask_new.nii.gz` in each subject's dwi/ folder.

# %% [Step 1a] Prepare group-specific atlases: N4-correct each atlas T2w
import time

# ============================================================
# GROUP-SPECIFIC ATLAS DEFINITIONS — one manually corrected subject per group
# ============================================================
GROUP_ATLAS_IDS = {
    "CTRL": "sub-129-921-R1F",
    "AP":   "sub-596-C-OO-NmM",
    "APTR": "sub-AP72-862-901-NmF",
}

# Build atlas info dict: {group: {t2w_n4, mask_nii, atlas_id, atlas_sub}}
atlas_info = {}
for group, atlas_id in GROUP_ATLAS_IDS.items():
    asub = next(s for s in subjects if s['id'] == atlas_id)
    info = {
        'atlas_id':   atlas_id,
        'atlas_sub':  asub,
        't2w':        asub['t2w'],
        'mask_mif':   asub['mask_new_edited'],
        'mask_nii':   asub['dwi_dir'] / "mask_new_edited.nii.gz",
        't2w_nii':    asub['dwi_dir'] / "_atlas_t2w.nii.gz",
        't2w_n4':     asub['dwi_dir'] / "_atlas_t2w_n4.nii.gz",
    }
    atlas_info[group] = info

print(f"{'=' * 70}")
print(f" GROUP-SPECIFIC ATLASES")
print(f"{'=' * 70}")
for group, info in atlas_info.items():
    print(f"  {group:>4}: {info['atlas_id']}")
    assert info['t2w'].exists(),     f"Atlas T2w not found: {info['t2w']}"
    assert info['mask_mif'].exists(), f"Atlas mask not found: {info['mask_mif']}"
    print(f"         T2w: ✅  Mask: ✅")

# N4-correct each atlas T2w and convert mask to NIfTI
for group, info in atlas_info.items():
    print(f"\nPreparing {group} atlas ({info['atlas_id']})...")
    run(f"{MODULES} && mrconvert {info['t2w']} {info['t2w_nii']} -force -quiet")
    run(f"{MODULES} && N4BiasFieldCorrection -d 3 -i {info['t2w_nii']} -o {info['t2w_n4']} "
        f"-s 2 -b [50] -c [200x200x200,0.0]")
    run(f"{MODULES} && mrconvert {info['mask_mif']} {info['mask_nii']} -force -quiet")
    check_file(info['t2w_n4'], f"{group} atlas T2w N4")
    check_file(info['mask_nii'], f"{group} atlas mask NIfTI")

print("\n✅ All 3 group atlases prepared.")

# %% [Step 1b] Register group atlas → each target and warp mask (PARALLEL)
from concurrent.futures import ProcessPoolExecutor, as_completed

# --- Worker function (runs in a separate process) ---
def register_and_warp_mask(args):
    """
    Register atlas T2w → target T2w using ANTs SyN, then warp the atlas mask
    into target space. Returns (sid, status, message).
    """
    import subprocess, os
    from pathlib import Path

    (sid, group, target_t2w_mif, dwi_dir_str, mask_nii_path, mask_mif_path,
     atlas_t2w_n4, atlas_mask_nii, modules) = args

    def _run(cmd, timeout=1800):
        return subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)

    dwi = Path(dwi_dir_str)
    target_t2w = Path(target_t2w_mif)
    mask_nii = Path(mask_nii_path)
    mask_mif = Path(mask_mif_path)

    if not target_t2w.exists():
        return (sid, 'skipped', 'T2w not found')
    if mask_mif.exists():
        return (sid, 'skipped', 'mask already exists')

    # Temp files
    prefix     = dwi / f"_tmp_{sid}_"               # ANTs output prefix
    tgt_nii    = dwi / f"_tmp_{sid}_t2w.nii.gz"
    tgt_n4     = dwi / f"_tmp_{sid}_t2w_n4.nii.gz"
    warped_nii = dwi / f"_tmp_{sid}_mask_warped.nii.gz"
    mask_clean = dwi / f"_tmp_{sid}_mask_clean.mif"

    # ANTs outputs: {prefix}Warped.nii.gz, {prefix}0GenericAffine.mat,
    # {prefix}1Warp.nii.gz, {prefix}1InverseWarp.nii.gz
    warp_file   = Path(f"{prefix}1Warp.nii.gz")
    affine_file = Path(f"{prefix}0GenericAffine.mat")

    steps = [
        # 1) Convert target T2w to NIfTI
        (f"{modules} && mrconvert {target_t2w} {tgt_nii} -force -quiet",
         "convert T2w to NIfTI"),

        # 2) N4 bias-correct target T2w
        (f"{modules} && N4BiasFieldCorrection -d 3 -i {tgt_nii} -o {tgt_n4} "
         f"-s 2 -b [50] -c [200x200x200,0.0]",
         "N4 bias correction"),

        # 3) Register atlas (moving) → target (fixed) using SyN
        (f"{modules} && antsRegistrationSyN.sh -d 3 "
         f"-f {tgt_n4} -m {atlas_t2w_n4} -o {prefix} -t s -n 4",
         "ANTs SyN registration"),

        # 4) Apply the transforms to warp atlas mask → target space
        (f"{modules} && antsApplyTransforms -d 3 "
         f"-i {atlas_mask_nii} -r {tgt_n4} -o {warped_nii} "
         f"-t {warp_file} -t {affine_file} "
         f"-n NearestNeighbor",
         "warp atlas mask"),

        # 5) Convert warped mask to MIF
        (f"{modules} && mrconvert {warped_nii} {mask_mif} -force -quiet",
         "convert mask to MIF"),

        # 6) Light cleanup: keep only largest connected component
        (f"{modules} && maskfilter {mask_mif} clean {mask_clean} -largest -force -quiet "
         f"&& mv {mask_clean} {mask_mif}",
         "morphological cleanup"),

        # 7) Save NIfTI copy
        (f"{modules} && mrconvert {mask_mif} {mask_nii} -force -quiet",
         "save NIfTI copy"),
    ]

    for cmd, step_name in steps:
        r = _run(cmd)
        if r.returncode != 0:
            # Cleanup temp files
            import glob as _glob
            for tmp in _glob.glob(str(prefix) + "*"):
                try: os.remove(tmp)
                except: pass
            for tmp in [tgt_nii, tgt_n4, warped_nii, mask_clean]:
                try: Path(tmp).unlink(missing_ok=True)
                except: pass
            return (sid, 'failed', f"{step_name} failed: {r.stderr[:200]}")

    # Cleanup ALL temp files (ANTs produces many)
    import glob as _glob
    for tmp in _glob.glob(str(prefix) + "*"):
        try: os.remove(tmp)
        except: pass
    for tmp in [tgt_nii, tgt_n4, warped_nii, mask_clean]:
        try: Path(tmp).unlink(missing_ok=True)
        except: pass

    # Get brain fraction
    frac = ""
    st = _run(f"{modules} && mrstats {mask_mif}")
    if st.returncode == 0:
        lines = st.stdout.strip().split('\n')
        if len(lines) >= 2:
            frac = f", brain fraction={lines[1].split()[1]}"

    return (sid, 'success', f"done{frac}")

# --- For atlas subjects, just copy their edited mask ---
def copy_atlas_mask(args):
    """Atlas subject: copy mask_new_edited → mask_new (no registration needed)."""
    import subprocess, shutil
    from pathlib import Path

    sid, _, _, _, mask_nii_path, mask_mif_path, _, atlas_mask_nii, modules = args
    mask_mif = Path(mask_mif_path)
    mask_nii = Path(mask_nii_path)
    edited = mask_mif.parent / "mask_new_edited.mif"

    if mask_mif.exists():
        return (sid, 'skipped', 'mask already exists')
    if not edited.exists():
        return (sid, 'failed', 'mask_new_edited.mif not found')

    shutil.copy2(str(edited), str(mask_mif))
    # Also save NIfTI copy
    r = subprocess.run(
        f"{modules} && mrconvert {mask_mif} {mask_nii} -force -quiet",
        shell=True, capture_output=True, text=True, timeout=60
    )
    if r.returncode != 0:
        return (sid, 'failed', f"NIfTI conversion failed: {r.stderr[:150]}")

    return (sid, 'success', 'copied from mask_new_edited.mif (atlas)')

# --- Dispatch parallel jobs ---
# ANTs SyN uses ~4 threads per job (-n 4). 9 targets total (3 per group).
# lmc02: 40 cores → 40/4 = 10 concurrent jobs (we have 9 targets + 3 atlas copies)
N_WORKERS = 9
ATLAS_IDS = set(GROUP_ATLAS_IDS.values())

print(f"{'=' * 70}")
print(f" GENERATING BRAIN MASKS (group-specific atlas registration)")
print(f" FOR {len(subjects)} SUBJECTS ({len(ATLAS_IDS)} atlases + {len(subjects) - len(ATLAS_IDS)} targets)")
print(f"{'=' * 70}")
for group, atlas_id in GROUP_ATLAS_IDS.items():
    targets = [s['id'] for s in subjects if s['group'] == group and s['id'] != atlas_id]
    print(f"  {group:>4} atlas: {atlas_id} → targets: {targets}")
print(f"\n Using {N_WORKERS} parallel workers on lmc02 (40 cores, 188 GB RAM)")
print(f" Registration: ANTs SyN (rigid + affine + nonlinear), ~3-5 min/subject")
print(f"{'=' * 70}")

t0 = time.time()
results = {'success': [], 'failed': [], 'skipped': []}

with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
    futures = {}
    for s in subjects:
        group = s['group']
        sid = s['id']
        info = atlas_info[group]

        args = (
            sid, group, str(s['t2w']), str(s['dwi_dir']),
            str(s['mask_new']), str(s['mask_new_mif']),
            str(info['t2w_n4']), str(info['mask_nii']), MODULES
        )

        if sid in ATLAS_IDS:
            # Atlas subject: just copy its edited mask
            futures[executor.submit(copy_atlas_mask, args)] = sid
        else:
            # Non-atlas: register group atlas → target and warp mask
            futures[executor.submit(register_and_warp_mask, args)] = sid

    for future in as_completed(futures):
        sid = futures[future]
        try:
            sid_ret, status, msg = future.result(timeout=1800)
            results[status].append(sid_ret)
            icon = {"success": "✅", "failed": "❌", "skipped": "⏭️"}[status]
            print(f"  {icon} {sid_ret}: {msg}")
        except Exception as e:
            results['failed'].append(sid)
            print(f"  ❌ {sid}: exception — {e}")

elapsed = time.time() - t0
print(f"\n{'=' * 70}")
print(f" MASK GENERATION RESULTS (group-specific atlas) — {elapsed:.0f}s total")
print(f"   ✅ Success: {len(results['success'])}")
print(f"   ❌ Failed:  {len(results['failed'])}")
print(f"   ⏭️  Skipped: {len(results['skipped'])}")
if results['failed']:
    print(f"   Failed subjects: {results['failed']}")
print(f"{'=' * 70}")

# %% [markdown] Step 2: Visual QC — Compare Old vs New Masks
# ---
# ## Step 2: Visual QC — Inspect Every Mask
#
# This step opens mrview for **each subject** with the atlas-based mask overlaid
# on the T2w image. Scroll through **every slice** and note any problems.
#
# **What to look for:**
# - Mask covers the **entire brain** without gaps
# - Mask **excludes skull/background**
# - For AP/APTR: check the **graft-host border** — the automatic mask may
#   include/exclude the graft incorrectly
#
# After reviewing, note which subjects need manual correction in the next step.

# %% [Step 2a] QC: View ONE subject to test (pick first CTRL)
# Quick test with the first subject to verify mrview works
#
# Copy-paste this into a BRODMANN126 terminal (absolute paths, no /export):
# mrview /home/glerma/glerma/PROYECTOS/mbrain/data/derivatives/mrtrix_newmask/sub-129-921-R1F/ses-260115121231/dwi/t2w_highres_flip.mif -overlay.load /home/glerma/glerma/PROYECTOS/mbrain/data/derivatives/mrtrix_newmask/sub-129-921-R1F/ses-260115121231/dwi/mask_new.mif -overlay.opacity 0.3 -overlay.colour 0,1,0
#
s = subjects[0]
print(f"Opening QC for: {s['id']} ({s['group']})")
if s['mask_new_mif'].exists():
    launch_gui(
        f"mrview {s['t2w']} "
        f"-overlay.load {s['mask_new_mif']} -overlay.opacity 0.3 "
        f"-overlay.colour 0,1,0"
    )
    print("Green overlay = atlas-based mask. Scroll through all slices.")
else:
    print(f"❌ Mask not found: {s['mask_new_mif']}")

# %% [Step 2b] QC: Compare OLD vs NEW mask side-by-side (first subject)
s = subjects[0]
print(f"Comparing OLD (red) vs NEW (green) mask for: {s['id']}")
if s['mask_new_mif'].exists() and s['mask_old'].exists():
    launch_gui(
        f"mrview {s['t2w']} "
        f"-overlay.load {s['mask_old']} -overlay.opacity 0.3 -overlay.colour 1,0,0 "
        f"-overlay.load {s['mask_new_mif']} -overlay.opacity 0.3 -overlay.colour 0,1,0"
    )
    print("RED = old dwi2mask | GREEN = new atlas-based")
else:
    print("❌ One or both masks missing")

# %% [Step 2c] QC: Systematic review — ALL subjects, one by one
# Run this cell, review one subject, then press Enter to continue to the next.
# Record any subjects that need manual editing below.

SUBJECTS_TO_FIX = []  # ← Will be populated as you review

for i, s in enumerate(subjects):
    label = f"[{i+1}/{len(subjects)}] {s['id']} ({s['group']})"
    
    if not s['mask_new_mif'].exists():
        print(f"  ⏭️ {label}: no mask found, skipping")
        continue
    
    print(f"\n{'=' * 60}")
    print(f"  QC: {label}")
    print(f"{'=' * 60}")
    
    launch_gui(
        f"mrview {s['t2w']} "
        f"-overlay.load {s['mask_new_mif']} -overlay.opacity 0.3 "
        f"-overlay.colour 0,1,0"
    )
    
    answer = input(f"  {s['id']}: PASS (Enter) or FAIL (type 'f'): ").strip().lower()
    if answer == 'f':
        SUBJECTS_TO_FIX.append(s['id'])
        print(f"  → Marked for manual correction: {s['id']}")
    else:
        print(f"  → PASS ✅")

print(f"\n{'=' * 60}")
print(f" QC COMPLETE")
print(f"   Subjects needing manual correction: {len(SUBJECTS_TO_FIX)}")
for sid in SUBJECTS_TO_FIX:
    print(f"     - {sid}")
print(f"{'=' * 60}")

# %% [markdown] Step 3: Manual Mask Correction
# ---
# ## Step 3: Manual Mask Correction (for failed QC subjects)
#
# For subjects where the automatic mask got the graft-host border wrong:
# 1. Open the mask in mrview's ROI editor
# 2. Paint/erase voxels to fix the boundary
# 3. Save as `mask_new_edited.mif`
#
# **mrview ROI editor usage:**
# - Tool → ROI editor
# - Load the atlas-based mask as the ROI
# - Left-click = add voxels, Right-click = remove voxels
# - Adjust brush size with scroll wheel
# - File → Save to save the edited mask

# %% [Step 3a] Manual Correction — Set subjects to fix
# ============================================================
# EDIT THIS: List subjects that failed QC in Step 2c.
# If Step 2c was run interactively, SUBJECTS_TO_FIX is already set.
# Otherwise, manually specify them here:
# ============================================================

# Uncomment and edit if needed:
# SUBJECTS_TO_FIX = ["sub-AP72-862-901-NmF", "sub-Apallial17NmF"]

print(f"Subjects to manually correct: {SUBJECTS_TO_FIX}")
if not SUBJECTS_TO_FIX:
    print("None! All masks passed QC. Skip to Step 4.")

# %% [Step 3b] Open mrview ROI editor for each subject needing correction
for sid in SUBJECTS_TO_FIX:
    s = next((x for x in subjects if x['id'] == sid), None)
    if s is None:
        print(f"❌ Subject {sid} not found in subject list!")
        continue
    
    src = s['mask_new_mif']
    dst = s['mask_new_edited']
    
    if not src.exists():
        print(f"❌ {sid}: mask not found")
        continue
    
    # Copy mask to edited version (so we don't lose the original)
    if not dst.exists():
        run(f"cp {src} {dst}", verbose=False)
    
    print(f"\n{'=' * 60}")
    print(f"  EDITING: {sid} ({s['group']})")
    print(f"  Mask file: {dst}")
    print(f"{'=' * 60}")
    print("  Instructions:")
    print("    1. In mrview: Tool → ROI editor")
    print("    2. Load ROI: the mask_new_edited.mif")
    print("    3. Paint (left-click) / Erase (right-click) to fix boundaries")
    print("    4. Save (Ctrl+S or File → Save) when done")
    print("    5. Close mrview and press Enter to continue")
    
    launch_gui(
        f"mrview {s['t2w']} "
        f"-overlay.load {dst} -overlay.opacity 0.4 -overlay.colour 0,1,0"
    )
    
    input(f"  Press Enter when done editing {sid}...")
    
    if dst.exists():
        print(f"  ✅ Saved: {dst}")
    else:
        print(f"  ⚠️ Edited mask not found at {dst}")

print("\n✅ Manual correction complete.")

# %% [markdown] Step 4: Apply New Masks — Replace Old Masks
# ---
# ## Step 4: Apply New Masks
#
# For each subject:
# - Back up old `mask.mif` and `mask_anat.mif`
# - **Resample** the new mask from T2w resolution (46×134×160, 0.3mm) to DWI
#   resolution (28×100×107, 0.5mm) using `mrgrid regrid -interp nearest`
# - Install the resampled mask as both `mask.mif` and `mask_anat.mif`
#
# **Why resample?** The new mask was drawn on T2w_highres_flip (0.3×0.1×0.1mm),
# but all downstream processing (FOD, tractography, fixels) operates at DWI
# resolution (0.5×0.15×0.15mm). MRtrix commands like dwi2fod, tckgen, and
# mtnormalise expect the mask to match the DWI grid exactly.
#
# **Why NOT re-run DWI preprocessing (step 4)?** Running batch_process_flip1.py
# --steps 4 would overwrite our manually edited mask with a fresh dwi2mask
# output. The DWI denoising, bias correction, and registration are already
# fine from the original pipeline — only the mask needs replacing.

# %% [Step 4] Back up old masks and install new ones (resampled to DWI space)
BACKUP_SUFFIX = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

for i, s in enumerate(subjects):
    label = f"[{i+1}/{len(subjects)}] {s['id']} ({s['group']})"
    
    # Choose best available mask: edited > atlas-based > skip
    if s['mask_new_edited'].exists():
        new_mask = s['mask_new_edited']
        source = "edited"
    elif s['mask_new_mif'].exists():
        new_mask = s['mask_new_mif']
        source = "atlas-based"
    else:
        print(f"  ⏭️ {label}: no new mask available, skipping")
        continue
    
    dwi = s['dwi_dir']
    mask_mif = dwi / "mask.mif"
    mask_anat = dwi / "mask_anat.mif"
    
    # Use the registered DWI as the template for resampling
    dwi_template = dwi / "dwi_flip_gradcheck_denoise_bias_anat.mif"
    if not dwi_template.exists():
        # Fallback: use existing mask or FOD as template (same grid)
        for fallback in ["wm_fod_norm.mif", "dwi_preproc.mif"]:
            fb = dwi / fallback
            if fb.exists():
                dwi_template = fb
                break
    
    if not dwi_template.exists():
        print(f"  ❌ {label}: no DWI template found for resampling, skipping")
        continue
    
    # Backup old masks
    if mask_mif.exists():
        backup = dwi / f"mask_old_{BACKUP_SUFFIX}.mif"
        run(f"mv {mask_mif} {backup}", verbose=False)
    if mask_anat.exists():
        backup = dwi / f"mask_anat_old_{BACKUP_SUFFIX}.mif"
        run(f"mv {mask_anat} {backup}", verbose=False)
    
    # Resample new mask from T2w resolution to DWI resolution
    mask_resampled = dwi / "mask_new_dwi_space.mif"
    result = run(
        f"{MODULES} && mrgrid {new_mask} regrid {mask_resampled} "
        f"-template {dwi_template} -interp nearest -force -quiet",
        verbose=False
    )
    
    if result.returncode != 0:
        print(f"  ❌ {label}: mrgrid resample FAILED")
        continue
    
    # Install resampled mask as both mask.mif and mask_anat.mif
    run(f"cp {mask_resampled} {mask_mif}", verbose=False)
    run(f"cp {mask_resampled} {mask_anat}", verbose=False)
    
    print(f"  ✅ {label}: installed {source} mask (resampled to DWI grid)")

print(f"\n✅ All masks updated. Old masks backed up with suffix: {BACKUP_SUFFIX}")

# %% [Step 4-verify] Verify installed masks match DWI grid
print(f"{'=' * 70}")
print(f" VERIFYING INSTALLED MASKS MATCH DWI GRID")
print(f"{'=' * 70}")

for i, s in enumerate(subjects):
    dwi = s['dwi_dir']
    mask_f = dwi / "mask.mif"
    fod_f = dwi / "wm_fod_norm.mif"
    
    if not mask_f.exists():
        print(f"  ❌ {s['id']}: mask.mif not found!")
        continue
    
    # Get dims via mrinfo
    r_mask = subprocess.run(f"mrinfo {mask_f} -size", shell=True, capture_output=True, text=True)
    mask_dims = r_mask.stdout.strip()
    
    if fod_f.exists():
        r_fod = subprocess.run(f"mrinfo {fod_f} -size", shell=True, capture_output=True, text=True)
        fod_dims = ' '.join(r_fod.stdout.strip().split()[:3])  # first 3 dims only
        match = "✅" if mask_dims == fod_dims else "❌ MISMATCH"
        print(f"  {match} {s['id']}: mask={mask_dims}  FOD={fod_dims}")
    else:
        print(f"  ⚠️  {s['id']}: mask={mask_dims}  (no FOD to compare)")

# %% [Step 4-checkpoint] Check pipeline readiness — what needs re-running?
# Run this BEFORE any Step 5 cell to see what still needs processing.
print(f"{'=' * 70}")
print(f" PIPELINE CHECKPOINT — WHAT NEEDS RE-RUNNING?")
print(f"{'=' * 70}")

chk = {'mask_ok': [], 'mask_bad': [],
       'fod_ok': [], 'fod_stale': [], 'fod_missing': [],
       'grp_ok': [], 'grp_stale': [], 'grp_missing': [],
       'tck_ok': [], 'tck_stale': [], 'tck_missing': []}

for s in subjects:
    dwi = s['dwi_dir']; sid = s['id']
    mask_f   = dwi / "mask.mif"
    mask_dws = dwi / "mask_new_dwi_space.mif"
    fod_f    = dwi / "wm_fod_norm.mif"
    grp_f    = dwi / "wm_fod_norm_group.mif"
    tck_f    = s['wb_tck']

    # 1. Mask at DWI resolution?
    if mask_f.exists():
        r = subprocess.run(f"mrinfo {mask_f} -size", shell=True, capture_output=True, text=True)
        dims = r.stdout.strip().split()
        chk['mask_ok' if (len(dims) == 3 and dims != ['46','134','160']) else 'mask_bad'].append(sid)
    else:
        chk['mask_bad'].append(sid)

    # 2. FOD newer than resampled mask?
    if fod_f.exists() and mask_dws.exists():
        chk['fod_ok' if fod_f.stat().st_mtime > mask_dws.stat().st_mtime else 'fod_stale'].append(sid)
    elif fod_f.exists():
        chk['fod_stale'].append(sid)
    else:
        chk['fod_missing'].append(sid)

    # 3. Group FOD newer than per-subject FOD?
    if grp_f.exists() and fod_f.exists():
        chk['grp_ok' if grp_f.stat().st_mtime > fod_f.stat().st_mtime else 'grp_stale'].append(sid)
    elif grp_f.exists():
        chk['grp_stale'].append(sid)
    else:
        chk['grp_missing'].append(sid)

    # 4. Tractogram newer than FOD?
    if tck_f.exists() and fod_f.exists():
        chk['tck_ok' if tck_f.stat().st_mtime > fod_f.stat().st_mtime else 'tck_stale'].append(sid)
    elif tck_f.exists():
        chk['tck_stale'].append(sid)
    else:
        chk['tck_missing'].append(sid)

print(f"\n  MASKS (DWI resolution):  ✅ {len(chk['mask_ok'])}/12")
if chk['mask_bad']:
    print(f"    ❌ BAD: {chk['mask_bad']}  → Run Step 4 first!")

needs_fod = chk['fod_stale'] + chk['fod_missing']
print(f"\n  FODs (wm_fod_norm.mif):  ✅ {len(chk['fod_ok'])}/12  stale/missing={len(needs_fod)}")
if needs_fod:
    print(f"    → Subjects needing FOD: {needs_fod}")

needs_grp = chk['grp_stale'] + chk['grp_missing']
print(f"\n  GROUP FODs:              ✅ {len(chk['grp_ok'])}/12  stale/missing={len(needs_grp)}")

needs_tck = chk['tck_stale'] + chk['tck_missing']
print(f"\n  TRACTOGRAMS:             ✅ {len(chk['tck_ok'])}/12  stale/missing={len(needs_tck)}")

print(f"\n{'=' * 70}")
print(f" RECOMMENDED NEXT STEP:")
if chk['mask_bad']:
    print(f"   → Step 4: Install masks (some missing/wrong resolution)")
elif needs_fod:
    print(f"   → Step 5b: FOD estimation ({len(needs_fod)} subjects)")
elif needs_grp:
    print(f"   → Step 5b-group: Group response & normalization")
elif needs_tck:
    print(f"   → Step 5d: Tractography ({len(needs_tck)} subjects)")
else:
    print(f"   → Everything up to date! Check Steps 5e+ for downstream.")
print(f"{'=' * 70}")

# %% [markdown] Step 5: Re-run Downstream Pipeline
# ---
# ## Step 5: Re-run Downstream Pipeline
#
# With new masks installed (resampled to DWI resolution), we re-run everything
# that depends on the mask. We skip DWI preprocessing (step 4 in
# batch_process_flip1.py) because:
# - It would overwrite our mask with a fresh dwi2mask output
# - The DWI denoising, bias correction, and registration are unchanged
# - The N4 bias correction quality difference from using a slightly different
#   mask boundary is negligible
#
# **Parallelization:** All per-subject steps use `run_parallel()` with up to
# 12 workers (one per subject) on lmc02's 40 cores.
#
# | Step | What | Why re-run |
# |------|------|------------|
# | 5b | FOD estimation | Uses mask for response estimation + CSD |
# | 5b-group | Group response + normalization | Uses mask for mtnormalise |
# | 5c | ROI resampling | FOD space reference |
# | 5d | Whole-brain tractography | Uses mask as seed/constraint |
# | 5e | SIFT2 | Depends on tractogram |
# | 5f | Passthrough extraction | Depends on tractogram |
# | 5g | MATLAB cleaning | Depends on tracts |
# | 5h | Metrics report | Depends on cleaned tracts |
# | 5i | Fixel analysis | Uses group-normalized FOD + mask |

# %% [Step 5b] Re-run FOD estimation (step 5) — PARALLEL (12 subjects)
print(f"{'=' * 70}")
print(f" RE-RUNNING FOD ESTIMATION WITH NEW MASKS (PARALLEL)")
print(f" Workers: {MAX_WORKERS} on {multiprocessing.cpu_count()} cores")
print(f"{'=' * 70}")

fod_tasks = [
    (s['id'],
     f"{MODULES} && cd {REPO_DIR} && python functions/batch_process_flip1.py "
     f"--subject {s['id']} --steps 5 --force --derivatives-dir {DERIV_BASE}",
     1800)
    for s in subjects
]
fod_results = run_parallel(fod_tasks, label="FOD estimation")

# %% [markdown] Step 5b-group: Group Response & Normalized FOD
# ---
# ## Step 5b-group: Recompute Group Response Functions & Group-Normalized FODs
#
# **Why this step is critical:**
# The fixel analysis (Step 5i) and any cross-subject AFD comparison use
# `wm_fod_norm_group.mif` — FODs estimated with a **group-averaged** response
# function, then normalized with `mtnormalise`. This ensures that AFD/FD values
# are comparable between subjects and between groups (CTRL vs AP vs APTR).
#
# Changing masks affects ALL of:
# 1. `dwi2response dhollander` — the mask constrains which voxels are used to
#    estimate per-subject response functions (WM, GM, CSF)
# 2. `dwi2fod msmt_csd` — the mask constrains where FODs are estimated
# 3. `mtnormalise` — the mask defines the brain volume for intensity normalization
#
# **Pipeline:**
# 1. Collect updated per-subject response functions (from Step 5b)
# 2. Average them across ALL subjects → `group_responses/group_response_{wm,gm,csf}.txt`
# 3. Re-estimate FODs per subject using the group response → `wm_fod_group.mif`
# 4. Re-normalize per subject → `wm_fod_norm_group.mif`
#
# If `mtnormalise` fails (known issue for some mouse data), the script falls
# back to using the unnormalized group FOD.

# %% [Step 5b-group] Recompute Group Response & Normalized FOD
print(f"{'=' * 70}")
print(f" RECOMPUTING GROUP RESPONSE FUNCTIONS & GROUP-NORMALIZED FODS")
print(f" (This uses ALL {len(subjects)} subjects for response averaging)")
print(f"{'=' * 70}")

result_grp = run(
    f"{MODULES} && cd {REPO_DIR} && python functions/compute_group_response.py "
    f"--force --derivatives-dir {DERIV_BASE}",
    timeout=7200
)

if result_grp.returncode == 0:
    print(f"\n  ✅ Group response & FOD normalization complete")
else:
    print(f"\n  ❌ Group response computation FAILED — check output above")

# Verify wm_fod_norm_group.mif exists for each subject
grp_fod_ok = 0
grp_fod_missing = 0
for s in subjects:
    fod_grp = s['dwi_dir'] / 'wm_fod_norm_group.mif'
    if fod_grp.exists():
        grp_fod_ok += 1
    else:
        grp_fod_missing += 1
        print(f"  ❌ MISSING: {s['id']} wm_fod_norm_group.mif")

print(f"\n  wm_fod_norm_group.mif: ✅ {grp_fod_ok} present, ❌ {grp_fod_missing} missing")

# %% [Step 5b-verify] Verify FOD estimation — dates, sizes, grids, normalisation
import datetime, json as _json
print(f"{'=' * 70}")
print(f" STEP 5b VERIFICATION — FOD Estimation & Group Normalisation")
print(f"{'=' * 70}\n")

_fod_files = [
    'wm_response.txt', 'gm_response.txt', 'csf_response.txt',  # per-subject responses
    'wm_fod.mif', 'gm_fod.mif', 'csf_fod.mif',                # per-subject FODs
    'wm_fod_norm.mif',                                          # per-subject normalised FOD
    'wm_fod_group.mif', 'gm_fod_group.mif', 'csf_fod_group.mif',  # group FODs
    'wm_fod_norm_group.mif',                                    # group-normalised FOD (key output)
]
_mask_files = ['mask.mif', 'mask_anat.mif', 'mask_new_dwi_space.mif']

_issues = []
_table_rows = []

for s in subjects:
    DWI = s['dwi_dir']
    sid = s['id']
    row = {'subject': sid, 'group': s['group']}

    # -- Check key output exists --
    fod_norm_grp = DWI / 'wm_fod_norm_group.mif'
    if not fod_norm_grp.exists():
        _issues.append(f"❌ {sid}: wm_fod_norm_group.mif MISSING")
        row['status'] = 'MISSING'
        _table_rows.append(row)
        continue

    # -- File modification time --
    mtime = datetime.datetime.fromtimestamp(fod_norm_grp.stat().st_mtime)
    row['mtime'] = mtime.strftime('%Y-%m-%d %H:%M')
    row['size_MB'] = round(fod_norm_grp.stat().st_size / 1e6, 1)

    # -- Check if normalised == unnormalised (mtnormalise fallback detection) --
    fod_grp_raw = DWI / 'wm_fod_group.mif'
    csf_fod_norm = DWI / 'csf_fod_norm_group.mif'
    if fod_grp_raw.exists() and fod_norm_grp.exists():
        sz_raw = fod_grp_raw.stat().st_size
        sz_norm = fod_norm_grp.stat().st_size
        # If sizes match exactly AND mtime within 1s → likely a raw copy (both fallbacks failed)
        if sz_raw == sz_norm and abs(fod_grp_raw.stat().st_mtime - fod_norm_grp.stat().st_mtime) < 2:
            row['normalised'] = '❌ RAW COPY (not normalised)'
            _issues.append(f"❌ {sid}: mtnormalise completely failed — wm_fod_norm_group.mif is a raw copy")
        elif not csf_fod_norm.exists() and sz_raw != sz_norm:
            row['normalised'] = '⚠️  WM-only norm'
            # WM-only is acceptable, just note it
        else:
            row['normalised'] = '✅ normalised'
    else:
        row['normalised'] = '? (cannot check)'

    # -- Grid check: FOD vs mask vs DWI --
    try:
        fod_sz = run(f"{MODULES} && mrinfo {fod_norm_grp} -size -quiet", timeout=30, verbose=False)
        mask_sz = run(f"{MODULES} && mrinfo {DWI / 'mask_anat.mif'} -size -quiet", timeout=30, verbose=False)
        fod_grid = fod_sz.stdout.strip()
        mask_grid = mask_sz.stdout.strip()
        row['fod_grid'] = fod_grid
        row['mask_grid'] = mask_grid
        # Compare first 3 dims only (FOD has 4th dim for SH coefficients)
        fod_dims = fod_grid.split()[:3]
        mask_dims = mask_grid.split()[:3]
        if fod_dims != mask_dims:
            _issues.append(f"❌ {sid}: GRID MISMATCH — FOD {fod_dims} vs mask {mask_dims}")
            row['grid_match'] = '❌ MISMATCH'
        else:
            row['grid_match'] = '✅'
    except Exception as e:
        row['grid_match'] = f'? error: {e}'

    # -- Check per-subject response files exist --
    missing_resp = [f for f in ['response_wm.txt', 'response_gm.txt', 'response_csf.txt']
                    if not (DWI / f).exists()]
    if missing_resp:
        _issues.append(f"⚠️  {sid}: missing response files: {missing_resp}")
        row['responses'] = f'⚠️  missing {len(missing_resp)}'
    else:
        row['responses'] = '✅'

    # -- Check mask_anat.mif is the new mask (should be recent) --
    mask_anat = DWI / 'mask_anat.mif'
    if mask_anat.exists():
        mask_mt = datetime.datetime.fromtimestamp(mask_anat.stat().st_mtime)
        row['mask_date'] = mask_mt.strftime('%Y-%m-%d %H:%M')
    else:
        row['mask_date'] = 'MISSING'
        _issues.append(f"❌ {sid}: mask_anat.mif MISSING")

    row['status'] = 'OK'
    _table_rows.append(row)

# -- Group response files --
grp_resp_dir = DERIV_BASE / 'group_responses'
print("  Group response files:")
for gf in ['group_response_wm.txt', 'group_response_gm.txt', 'group_response_csf.txt']:
    gp = grp_resp_dir / gf
    if gp.exists():
        gmt = datetime.datetime.fromtimestamp(gp.stat().st_mtime).strftime('%Y-%m-%d %H:%M')
        print(f"    ✅ {gf}  (modified: {gmt})")
    else:
        print(f"    ❌ {gf}  MISSING")
        _issues.append(f"❌ Group response {gf} MISSING")

# -- Print per-subject table --
print(f"\n  {'Subject':<28} {'Group':<6} {'Status':<8} {'Norm?':<30} {'FOD MB':<8} "
      f"{'FOD Grid':<20} {'Grid?':<4} {'Resp':<4} {'FOD Date':<16} {'Mask Date':<16}")
print(f"  {'─'*28} {'─'*6} {'─'*8} {'─'*30} {'─'*8} {'─'*20} {'─'*4} {'─'*4} {'─'*16} {'─'*16}")
for r in _table_rows:
    print(f"  {r.get('subject','?'):<28} {r.get('group','?'):<6} {r.get('status','?'):<8} "
          f"{r.get('normalised','?'):<30} {r.get('size_MB','?'):<8} "
          f"{r.get('fod_grid','?'):<20} {r.get('grid_match','?'):<4} "
          f"{r.get('responses','?'):<4} {r.get('mtime','?'):<16} {r.get('mask_date','?'):<16}")

# -- Summary --
print(f"\n  {'=' * 60}")
if _issues:
    print(f"  ⚠️  {len(_issues)} issue(s) found:")
    for iss in _issues:
        print(f"    {iss}")
else:
    print(f"  ✅ ALL {len(subjects)} subjects passed Step 5b verification")
print(f"  {'=' * 60}")

# %% [Step 5c] Re-run ROI resampling + segmentation (steps 6,7) — PARALLEL
print(f"{'=' * 70}")
print(f" RE-RUNNING ROI RESAMPLING + SEGMENTATION (PARALLEL)")
print(f"{'=' * 70}")

seg_tasks = [
    (s['id'],
     f"{MODULES} && cd {REPO_DIR} && python functions/batch_process_flip1.py "
     f"--subject {s['id']} --steps 6,7 --force --derivatives-dir {DERIV_BASE}",
     1800)
    for s in subjects
]
seg_results = run_parallel(seg_tasks, label="ROI resampling")

# %% [markdown] Step 5c-QC: ROI Space & Laterality Verification
# ---
# ## Step 5c-QC: ROI Space & Laterality Verification
#
# After ROI resampling we must verify:
# 1. **Correct space** — fod_space ROI grid matches the FOD (dims & voxel size)
# 2. **Non-empty** — every ROI has >0 voxels after resampling
# 3. **Left/Right separation** — left and right ROIs have centroids on opposite
#    sides of axis 2 (L/R axis after flip1), with no overlap
# 4. **Within-subject consistency** — all "left" ROIs are on the same side,
#    all "right" ROIs on the opposite side
#
# **Axis convention (flip1):**
# - axis 0 = AP, axis 1 = DV, axis 2 = LR
# - CTRL (atlas):  left → LOW axis-2 indices,  right → HIGH
# - AP/APTR (manual+flip): left → HIGH axis-2 indices, right → LOW

# %% [Step 5c-QC] ROI Space & Laterality Check
import subprocess as _sp, tempfile as _tf, json as _json
import numpy as np

def _load_roi_data(roi_path):
    """Load a .mif or .nii.gz ROI and return (data_array, header_dict)."""
    p = str(roi_path)
    if p.endswith('.mif'):
        tmp = _tf.mktemp(suffix='.nii.gz')
        _sp.run(['mrconvert', p, tmp, '-quiet', '-force'],
                capture_output=True, timeout=30)
        import nibabel as nib
        img = nib.load(tmp)
        data = img.get_fdata()
        os.unlink(tmp)
    else:
        import nibabel as nib
        img = nib.load(p)
        data = img.get_fdata()
    return data

def _get_grid_info(mif_path):
    """Return (dims_3d, voxel_sizes_3d) for a MRtrix/NIfTI image."""
    r = _sp.run(
        f"mrinfo {mif_path} -size -spacing",
        shell=True, capture_output=True, text=True, timeout=10
    )
    lines = r.stdout.strip().split('\n')
    dims = [int(x) for x in lines[0].split()][:3]  # first 3 spatial dims
    vox = [float(x) for x in lines[1].split()][:3]
    return dims, vox

print(f"{'=' * 80}")
print(f" ROI SPACE & LATERALITY VERIFICATION — {len(subjects)} SUBJECTS")
print(f"{'=' * 80}")

roi_qc_issues = []  # collect (subject, issue_str)

for si, s in enumerate(subjects):
    label = f"[{si+1}/{len(subjects)}] {s['id']} ({s['group']})"
    print(f"\n{'─' * 70}")
    print(f"  {label}")
    print(f"{'─' * 70}")

    DWI = s['dwi_dir']
    fod_file = DWI / 'wm_fod_norm.mif'
    mask_file = DWI / 'mask.mif'

    # --- 1. Reference grid from FOD ---
    if not fod_file.exists():
        msg = "FOD not found — skipping"
        print(f"  ❌ {msg}")
        roi_qc_issues.append((s['id'], msg))
        continue

    fod_dims, fod_vox = _get_grid_info(fod_file)
    print(f"  FOD grid:  {fod_dims}  voxel: {[f'{v:.4f}' for v in fod_vox]}")

    # --- 2. Enumerate fod_space ROIs ---
    if s['group'] == 'CTRL':
        roi_base = s['base'] / 'atlas' / 'rois' / 'fod_space'
        roi_pairs = {
            'isocortex': ('roi_isocortex_left.nii.gz', 'roi_isocortex_right.nii.gz'),
            'thalamus':  ('roi_thalamus_left.nii.gz',  'roi_thalamus_right.nii.gz'),
            'striatum':  ('roi_striatum_left.nii.gz',  'roi_striatum_right.nii.gz'),
            'cerebellum':('roi_cerebellum_left.nii.gz', 'roi_cerebellum_right.nii.gz'),
        }
    else:
        roi_base = s['base'] / 'manual_roi' / 'fod_space'
        roi_pairs = {
            'Graft':     ('roi_Graft_left_flip_fod.mif',     'roi_Graft_right_flip_fod.mif'),
            'Thalamus':  ('roi_Thalamus_left_flip_fod.mif',  'roi_Thalamus_right_flip_fod.mif'),
            'Striatum':  ('roi_Striatum_left_flip_fod.mif',  'roi_Striatum_right_flip_fod.mif'),
            'Palecortex':('roi_Palecortex_left_flip_fod.mif','roi_Palecortex_right_flip_fod.mif'),
            'cerebellum':('roi_cerebellum_left_flip_fod.mif','roi_cerebellum_right_flip_fod.mif'),
        }

    left_centroids_ax2 = []
    right_centroids_ax2 = []

    for region, (lf, rf) in roi_pairs.items():
        for hemi_tag, fname in [('L', lf), ('R', rf)]:
            roi_path = roi_base / fname
            if not roi_path.exists():
                msg = f"{fname}: NOT FOUND"
                print(f"    ❌ {msg}")
                roi_qc_issues.append((s['id'], msg))
                continue

            # Check grid match
            roi_dims, roi_vox = _get_grid_info(roi_path)
            grid_ok = (roi_dims == fod_dims)
            vox_ok = all(abs(rv - fv) < 0.01 for rv, fv in zip(roi_vox, fod_vox))

            if not grid_ok:
                msg = f"{fname}: GRID MISMATCH dims={roi_dims} vs FOD={fod_dims}"
                print(f"    ❌ {msg}")
                roi_qc_issues.append((s['id'], msg))
            if not vox_ok:
                msg = f"{fname}: VOXEL SIZE MISMATCH vox={roi_vox} vs FOD={fod_vox}"
                print(f"    ⚠️  {msg}")
                roi_qc_issues.append((s['id'], msg))

            # Load data and check non-empty + centroid
            data = _load_roi_data(roi_path)
            coords = np.argwhere(data > 0.5)
            nvox = len(coords)

            if nvox == 0:
                msg = f"{fname}: EMPTY (0 voxels)"
                print(f"    ❌ {msg}")
                roi_qc_issues.append((s['id'], msg))
                continue

            centroid = coords.mean(axis=0)
            mid2 = data.shape[2] / 2.0
            side_str = "HIGH" if centroid[2] > mid2 else "LOW"

            status = "✅" if grid_ok else "⚠️ "
            print(f"    {status} {region:>12}_{hemi_tag}: {nvox:>5} vox  "
                  f"centroid_ax2={centroid[2]:5.1f}  mid={mid2:.1f}  side={side_str}")

            if hemi_tag == 'L':
                left_centroids_ax2.append((region, centroid[2]))
            else:
                right_centroids_ax2.append((region, centroid[2]))

    # --- 3. Left/Right consistency check ---
    if left_centroids_ax2 and right_centroids_ax2:
        # All lefts should be on the same side, all rights on the opposite side
        left_sides  = ['HIGH' if c > fod_dims[2]/2.0 else 'LOW' for _, c in left_centroids_ax2]
        right_sides = ['HIGH' if c > fod_dims[2]/2.0 else 'LOW' for _, c in right_centroids_ax2]

        left_consistent = len(set(left_sides)) == 1
        right_consistent = len(set(right_sides)) == 1
        lr_separated = left_consistent and right_consistent and (left_sides[0] != right_sides[0])

        if lr_separated:
            print(f"\n    ✅ L/R laterality: all L={left_sides[0]}, all R={right_sides[0]} — CONSISTENT")
        else:
            if not left_consistent:
                msg = f"LEFT ROIs on mixed sides: {list(zip([r for r,_ in left_centroids_ax2], left_sides))}"
                print(f"\n    ❌ {msg}")
                roi_qc_issues.append((s['id'], msg))
            if not right_consistent:
                msg = f"RIGHT ROIs on mixed sides: {list(zip([r for r,_ in right_centroids_ax2], right_sides))}"
                print(f"\n    ❌ {msg}")
                roi_qc_issues.append((s['id'], msg))
            if left_consistent and right_consistent and left_sides[0] == right_sides[0]:
                msg = f"L and R centroids BOTH on {left_sides[0]} side — SWAP ERROR?"
                print(f"\n    ❌ {msg}")
                roi_qc_issues.append((s['id'], msg))

        # Check L/R overlap for each paired region
        for region, (lf, rf) in roi_pairs.items():
            lp = roi_base / lf
            rp = roi_base / rf
            if lp.exists() and rp.exists():
                ld = _load_roi_data(lp)
                rd = _load_roi_data(rp)
                overlap = np.sum((ld > 0.5) & (rd > 0.5))
                if overlap > 0:
                    msg = f"{region}: L/R OVERLAP = {overlap} voxels"
                    print(f"    ⚠️  {msg}")
                    roi_qc_issues.append((s['id'], msg))

# --- Summary ---
print(f"\n{'=' * 80}")
print(f" ROI QC SUMMARY")
print(f"{'=' * 80}")
if roi_qc_issues:
    print(f"  ⚠️  {len(roi_qc_issues)} issues found:")
    for subj, issue in roi_qc_issues:
        print(f"    {subj:<30} {issue}")
else:
    print(f"  ✅ All ROIs passed: correct space, non-empty, L/R consistent")
print(f"{'=' * 80}")

# %% [Step 5d] Re-run whole-brain tractography — PARALLEL
# NOTE: tckgen is CPU-intensive. With 12 subjects × ~3 threads each = 36 cores used.
# lmc02 has 40 cores so this saturates the machine — that's fine, it's the fastest way.
# Uses generate_wholebrain_tractogram.py (1M streamlines, iFOD2, mask-seeded)
print(f"{'=' * 70}")
print(f" RE-RUNNING WHOLE-BRAIN TRACTOGRAPHY (PARALLEL)")
print(f" 12 subjects × tckgen 1M streamlines on {multiprocessing.cpu_count()} cores")
print(f"{'=' * 70}")

tck_tasks = [
    (s['id'],
     f"{MODULES} && cd {REPO_DIR} && python functions/generate_wholebrain_tractogram.py "
     f"--subject {s['id']} --force --derivatives-dir {DERIV_BASE}",
     None)
    for s in subjects
]
tck_results = run_parallel(tck_tasks, label="Tractography")

# %% [Step 5d-post] Downsample to 10k + convert to .trk
# Downsample wholebrain_1000k.tck → wholebrain_10k.tck using tckedit -number
# Then convert wholebrain_10k.tck → wholebrain_10k.trk using nibabel
import numpy as np
import tempfile

print(f"{'=' * 70}")
print(f" DOWNSAMPLING TO 10k + CONVERTING TO .trk — {len(subjects)} SUBJECTS")
print(f"{'=' * 70}")

# --- Part 1: Downsample with tckedit (parallel) ---
ds_tasks = [
    (s['id'],
     f"{MODULES} && tckedit "
     f"{DERIV_BASE}/{s['id']}/{s['session']}/dwi/wholebrain/wholebrain_1000k.tck "
     f"{DERIV_BASE}/{s['id']}/{s['session']}/dwi/wholebrain/wholebrain_10k.tck "
     f"-number 10000 -force",
     None)
    for s in subjects
]
ds_results = run_parallel(ds_tasks, label="Downsample 10k")

# --- Part 2: Convert .tck → .trk using nibabel ---
import nibabel as nib
from nibabel.streamlines import load as tck_load
from nibabel.streamlines.trk import TrkFile

print(f"\n  Converting .tck → .trk ...")
trk_ok, trk_fail = 0, 0
for s in subjects:
    sid = s['id']
    ses = s['session']
    dwi_dir = DERIV_BASE / sid / ses / "dwi"
    wb_dir = dwi_dir / "wholebrain"
    tck_10k = wb_dir / "wholebrain_10k.tck"
    trk_10k = wb_dir / "wholebrain_10k.trk"
    t2w_mif = dwi_dir / "t2w_highres_flip.mif"

    if not tck_10k.exists():
        print(f"  ❌ {sid}: wholebrain_10k.tck not found")
        trk_fail += 1
        continue

    try:
        # Convert t2w .mif → temp .nii.gz for nibabel header
        with tempfile.NamedTemporaryFile(suffix=".nii.gz", delete=False) as tmp:
            tmp_nii = tmp.name
        subprocess.run(
            f"{MODULES} && mrconvert {t2w_mif} {tmp_nii} -force -quiet",
            shell=True, check=True, capture_output=True, timeout=30
        )
        ref_img = nib.load(tmp_nii)
        aff = ref_img.affine
        aff_inv = np.linalg.inv(aff)
        zooms = np.array(ref_img.header.get_zooms()[:3], dtype=np.float32)
        dims = np.array(ref_img.shape[:3], dtype=np.int16)

        # Load .tck streamlines
        tck = tck_load(str(tck_10k))
        streamlines = list(tck.streamlines)

        # World (RAS mm) → voxel coordinates
        vox_streamlines = []
        for sl in streamlines:
            ones = np.ones((sl.shape[0], 1))
            vox = (aff_inv @ np.hstack([sl, ones]).T).T[:, :3]
            vox_streamlines.append(vox.astype(np.float32))

        # Build TrkFile and save
        hdr = TrkFile.create_empty_header()
        hdr['dimensions'] = dims
        hdr['voxel_sizes'] = zooms
        hdr['voxel_order'] = b'RAS'
        hdr['vox_to_ras'] = aff.astype(np.float32)

        tractogram = nib.streamlines.Tractogram(vox_streamlines, affine_to_rasmm=aff)
        nib.streamlines.save(TrkFile(tractogram, header=hdr), str(trk_10k))
        trk_ok += 1
        print(f"  ✅ {sid}: {trk_10k.name}")

        # Cleanup temp file
        os.remove(tmp_nii)
    except Exception as e:
        trk_fail += 1
        print(f"  ❌ {sid}: {e}")

print(f"\n  .trk conversion: ✅ {trk_ok} / ❌ {trk_fail}")

# --- QC: visualize in mrview ---
# Run these on BRODMANN126 to visually check the 10k tractograms:
#
# D=/home/glerma/glerma/PROYECTOS/mbrain/data/derivatives/mrtrix_newmask
# for d in $D/sub-*/ses-*/dwi/wholebrain; do
#   sub=$(basename $(dirname $(dirname $(dirname "$d"))))
#   mrview $(dirname "$d")/t2w_highres_flip.mif \
#     -tractography.load "$d/wholebrain_10k.tck" &
#   sleep 0.5
# done

# %% [markdown] Step 5e: SIFT2 Weights
# ---
# ## Step 5e: SIFT2 Weights for Whole-Brain Tractogram
#
# Compute SIFT2 cross-section multipliers so each streamline gets a weight
# proportional to how well the tractogram matches the underlying FOD.

# %% [Step 5e] Compute SIFT2 Weights — PARALLEL
print(f"{'=' * 70}")
print(f" COMPUTING SIFT2 WEIGHTS (PARALLEL) — {len(subjects)} SUBJECTS")
print(f"{'=' * 70}")

# Check timestamps: if sift2_weights.txt is older than the tractogram, remove it
# so the pipeline script will recompute it instead of skipping.
import os
stale_count = 0
for s in subjects:
    wb_tck = s['wb_tck']
    wf = s['sift2_weights']
    if wf.exists() and wb_tck.exists():
        tck_mtime = os.path.getmtime(wb_tck)
        wf_mtime  = os.path.getmtime(wf)
        if wf_mtime < tck_mtime:
            print(f"  STALE: {s['id']} — weights ({datetime.datetime.fromtimestamp(wf_mtime):%Y-%m-%d %H:%M}) "
                  f"older than tractogram ({datetime.datetime.fromtimestamp(tck_mtime):%Y-%m-%d %H:%M}), removing")
            wf.unlink()
            mu = s.get('sift2_mu') or (wf.parent / "sift2_mu.txt")
            if mu.exists():
                mu.unlink()
            stale_count += 1
        else:
            print(f"  OK:    {s['id']} — weights ({datetime.datetime.fromtimestamp(wf_mtime):%Y-%m-%d %H:%M}) "
                  f"up-to-date vs tractogram ({datetime.datetime.fromtimestamp(tck_mtime):%Y-%m-%d %H:%M})")
    elif not wf.exists() and wb_tck.exists():
        print(f"  NEW:   {s['id']} — no weights file, will compute")
        stale_count += 1
    elif not wb_tck.exists():
        print(f"  SKIP:  {s['id']} — no tractogram found")
print(f"\n  → {stale_count} subject(s) need SIFT2 (re)computation\n")

sift2_tasks = [
    (s['id'],
     f"{MODULES} && cd {REPO_DIR} && python functions/run_passthrough_pipeline.py "
     f"--subject {s['id']} --step sift2 --derivatives-dir {DERIV_BASE}",
     1800)
    for s in subjects
]
sift2_results = run_parallel(sift2_tasks, label="SIFT2")

# %% [markdown] Step 5e-CST: Create & Edit CST ROIs for AP/APTR Subjects
# ---
# ## Step 5e-CST: CST (Corticospinal Tract) ROI Preparation
#
# Before passthrough extraction, AP and APTR subjects need **CST ROIs** (left & right)
# for the corticospinal tract analysis. These are manually drawn on the t2w_highres_flip
# image using mrview's ROI editor, saved to `manual_roi_flip/`, then regridded to FOD space.
#
# Source dir:  `PROYECTOS/mbrain/local/manual_roi_flip/`
# Target dir:  `<subject>/dwi/wholebrain/passthrough_cst/`
#
# Workflow:
# 1. Check which AP/APTR subjects already have CST ROIs
# 2. Create empty template ROIs for subjects that are missing them
# 3. Open mrview for manual drawing (on BROADMANN126)
# 4. Regrid finished ROIs to FOD space

# %% [Step 5e-CST-check] Check CST ROI status for AP/APTR subjects
import os
from pathlib import Path

# Fix MRtrix tmp dir (system default points to non-existent /flywheel/v0/output/tmp)
os.environ['MRTRIX_TMPDIR'] = '/tmp'

# CST ROIs live in each subject's manual_roi/ dir as roi_CST_{hemi}_flip.mif
# (same location & naming scheme as roi_Graft_left_flip.mif, roi_Thalamus_left_flip.mif, etc.)
# Legacy source dir (for the two subjects already drawn on Mac):
CST_LEGACY_DIR = Path("/home/glerma/glerma/PROYECTOS/mbrain/local/manual_roi_flip")

print(f"{'=' * 70}")
print(f" CST ROI STATUS CHECK — AP & APTR SUBJECTS")
print(f"{'=' * 70}")
print(f"  ROI naming: manual_roi/roi_CST_{{left,right}}_flip.mif")
print(f"  FOD-space:  wholebrain/passthrough_cst/roi_CST_{{left,right}}_fod.mif")
print()

cst_missing = []   # subjects that need CST ROIs created
cst_exists  = []   # subjects that already have CST ROIs in FOD space
cst_needs_regrid = []  # subjects that have source ROIs but no FOD-space ROIs

for s in subjects:
    if s['group'] == 'CTRL':
        continue  # CTRL subjects don't need CST ROIs

    sid = s['id']
    manual_roi_dir = s['base'] / 'manual_roi'
    cst_dir = s['wb_dir'] / "passthrough_cst"

    # Check FOD-space ROIs (final product)
    fod_left  = cst_dir / "roi_CST_left_fod.mif"
    fod_right = cst_dir / "roi_CST_right_fod.mif"
    has_fod = fod_left.exists() and fod_right.exists()

    # Check source ROIs in subject's manual_roi/ (new standard location)
    src_left  = manual_roi_dir / "roi_CST_left_flip.mif"
    src_right = manual_roi_dir / "roi_CST_right_flip.mif"
    has_src = src_left.exists() and src_right.exists()

    # Also check legacy location (manual_roi_flip/ shared dir) if not in manual_roi/
    if not has_src:
        for pat_l, pat_r in [
            (CST_LEGACY_DIR / f"roi_CST_left_{sid}.mif",     CST_LEGACY_DIR / f"roi_CST_right_{sid}.mif"),
            (CST_LEGACY_DIR / f"roi_CST_left_{sid}_CST.mif", CST_LEGACY_DIR / f"roi_CST_right_{sid}_CST.mif"),
        ]:
            if pat_l.exists(): src_left = pat_l
            if pat_r.exists(): src_right = pat_r
        has_src = src_left.exists() and src_right.exists()

    status = ""
    if has_fod:
        status = "✅ FOD-space ROIs exist"
        cst_exists.append(s)
    elif has_src:
        status = f"⚠️  Source ROIs exist ({src_left.parent.name}/), need regrid to FOD space"
        cst_needs_regrid.append((s, src_left, src_right))
    else:
        status = "❌ No CST ROIs — need to create & draw"
        cst_missing.append(s)

    print(f"  [{s['group']:>4}] {sid:<30} {status}")

print(f"\n  Summary:")
print(f"    ✅ Complete (FOD-space):  {len(cst_exists)}")
print(f"    ⚠️  Need regrid only:     {len(cst_needs_regrid)}")
print(f"    ❌ Need creation & draw:  {len(cst_missing)}")

# %% [Step 5e-CST-create] Create empty CST ROI templates for subjects that need them
# Creates empty bitwise ROIs in t2w_highres_flip space in each subject's manual_roi/
# for manual drawing in mrview. Naming: roi_CST_{hemi}_flip.mif
os.environ['MRTRIX_TMPDIR'] = '/tmp'

print(f"{'=' * 70}")
print(f" CREATING EMPTY CST ROI TEMPLATES FOR MANUAL DRAWING")
print(f"{'=' * 70}")

created_files = []
for s in cst_missing:
    sid = s['id']
    t2w = s['t2w']
    manual_roi_dir = s['base'] / 'manual_roi'

    if not t2w.exists():
        print(f"  ❌ {sid}: t2w not found at {t2w}")
        continue

    manual_roi_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n  {sid} ({s['group']}):")

    for hemi in ['left', 'right']:
        out_file = manual_roi_dir / f"roi_CST_{hemi}_flip.mif"
        if out_file.exists():
            print(f"    SKIP: {out_file.name} already exists")
            continue
        # Create an empty (all-zero) bitwise ROI matching t2w grid
        cmd = f"mrcalc {t2w} 0 -mult {out_file} -datatype bit -force -quiet"
        print(f"    Creating: {out_file.name} ... ", end="", flush=True)
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
        if result.returncode == 0 and out_file.exists():
            created_files.append(out_file)
            print(f"✅ ({out_file.stat().st_size} bytes)")
        else:
            print(f"❌ {result.stderr[:200]}")

print(f"\n  Created {len(created_files)} empty ROI templates")
print(f"  Now open them in mrview (next cell) to draw the CST ROIs manually.")

# %% [Step 5e-CST-edit] Open mrview for manual CST ROI drawing — run on BROADMANN126
# Opens mrview with:
#   - t2w_highres_flip.mif as the base image
#   - Existing manual ROIs as overlays for reference (Graft, Thalamus)
#   - Empty CST ROI loaded in the ROI editor for drawing
#
# Instructions:
#   1. Tool → ROI editor
#   2. The CST ROI is pre-loaded — select it in the ROI editor panel
#   3. Draw the CST region (left or right) on coronal/axial slices
#   4. SAVE (Ctrl+S or File → Save) before closing
#   5. Repeat for the other hemisphere
print(f"{'=' * 70}")
print(f" MANUAL CST ROI DRAWING — mrview on BROADMANN126")
print(f"{'=' * 70}")
print()
print("INSTRUCTIONS:")
print("  1. In mrview: Tool → ROI editor")
print("  2. The CST ROI is pre-loaded in the editor")
print("  3. Draw the CST region on coronal/axial slices")
print("  4. SAVE before closing mrview")
print("  5. Repeat for each hemisphere & subject")
print()

# Collect all AP/APTR subjects needing manual drawing
all_need_draw = []
for s in subjects:
    if s['group'] == 'CTRL':
        continue
    manual_roi_dir = s['base'] / 'manual_roi'
    cst_dir = s['wb_dir'] / "passthrough_cst"
    # Skip if FOD-space ROIs already exist
    if (cst_dir / "roi_CST_left_fod.mif").exists() and (cst_dir / "roi_CST_right_fod.mif").exists():
        continue
    # Check if source ROIs exist (either in manual_roi/ or legacy dir)
    for hemi in ['left', 'right']:
        roi_file = manual_roi_dir / f"roi_CST_{hemi}_flip.mif"
        if not roi_file.exists():
            # Check legacy location
            for legacy in [CST_LEGACY_DIR / f"roi_CST_{hemi}_{s['id']}.mif",
                           CST_LEGACY_DIR / f"roi_CST_{hemi}_{s['id']}_CST.mif"]:
                if legacy.exists():
                    roi_file = legacy
                    break
        if roi_file.exists():
            all_need_draw.append((s, hemi, roi_file))

for s, hemi, roi_file in all_need_draw:
    sid = s['id']
    t2w = s['t2w']

    # Build overlay list: show existing manual ROIs for reference
    overlays = ""
    manual_roi_dir = s['base'] / 'manual_roi'
    for roi_name in ['roi_Graft_left_flip.mif', 'roi_Graft_right_flip.mif',
                     'roi_Thalamus_left_flip.mif', 'roi_Thalamus_right_flip.mif']:
        roi_path = manual_roi_dir / roi_name
        if roi_path.exists():
            overlays += f" -overlay.load {roi_path} -overlay.opacity 0.3"

    print(f"\n{'=' * 60}")
    print(f"  EDITING: {sid} ({s['group']}) — {hemi} CST")
    print(f"  ROI file: {roi_file}")
    print(f"{'=' * 60}")
    print("  Instructions:")
    print("    1. In mrview: Tool → ROI editor")
    print("    2. The CST ROI is pre-loaded — select it in the editor panel")
    print("    3. Draw the CST region on coronal/axial slices")
    print("    4. Save (Ctrl+S or File → Save) when done")
    print("    5. Close mrview and press Enter to continue")

    launch_gui(f"mrview {t2w}{overlays} -roi.load {roi_file} -roi.opacity 0.5")

    input(f"  Press Enter when done editing {sid} {hemi} CST...")

    if roi_file.exists() and roi_file.stat().st_size > 1000:
        print(f"  ✅ Saved: {roi_file.name} ({roi_file.stat().st_size} bytes)")
    else:
        print(f"  ⚠️  ROI file may still be empty — check before regridding")

if not all_need_draw:
    print("  All AP/APTR subjects already have CST ROIs — nothing to draw!")
else:
    print(f"\n  ✅ Finished launching mrview for {len(all_need_draw)} CST ROIs.")
    print("  Run the next cell to regrid them to FOD space.")

# %% [Step 5e-CST-regrid] Regrid CST ROIs to FOD space and copy to passthrough_cst/
# After manual drawing, regrid source ROIs from t2w space to FOD space using
# mrgrid regrid with the subject's wm_fod_norm_group.mif as template.
os.environ['MRTRIX_TMPDIR'] = '/tmp'

print(f"{'=' * 70}")
print(f" REGRIDDING CST ROIs TO FOD SPACE")
print(f"{'=' * 70}")

all_to_regrid = []
for s in subjects:
    if s['group'] == 'CTRL':
        continue
    sid = s['id']
    cst_dir = s['wb_dir'] / "passthrough_cst"
    fod_left  = cst_dir / "roi_CST_left_fod.mif"
    fod_right = cst_dir / "roi_CST_right_fod.mif"
    if fod_left.exists() and fod_right.exists():
        print(f"  SKIP: {sid} — FOD-space CST ROIs already exist")
        continue
    # Find source ROIs: check manual_roi/ first, then legacy dir
    manual_roi_dir = s['base'] / 'manual_roi'
    src_left = src_right = None
    # Standard location
    sl = manual_roi_dir / "roi_CST_left_flip.mif"
    sr = manual_roi_dir / "roi_CST_right_flip.mif"
    if sl.exists(): src_left = sl
    if sr.exists(): src_right = sr
    # Legacy location fallback
    if not src_left or not src_right:
        for pat_l, pat_r in [
            (CST_LEGACY_DIR / f"roi_CST_left_{sid}.mif",     CST_LEGACY_DIR / f"roi_CST_right_{sid}.mif"),
            (CST_LEGACY_DIR / f"roi_CST_left_{sid}_CST.mif", CST_LEGACY_DIR / f"roi_CST_right_{sid}_CST.mif"),
        ]:
            if not src_left and pat_l.exists(): src_left = pat_l
            if not src_right and pat_r.exists(): src_right = pat_r

    if not src_left or not src_right:
        print(f"  ❌ {sid}: Source ROIs not found (checked manual_roi/ and legacy dir)")
        continue
    all_to_regrid.append((s, src_left, src_right))

regrid_ok = 0
for s, src_left, src_right in all_to_regrid:
    sid = s['id']
    fod_template = s['dwi_dir'] / "wm_fod_norm_group.mif"
    cst_dir = s['wb_dir'] / "passthrough_cst"
    cst_dir.mkdir(parents=True, exist_ok=True)

    if not fod_template.exists():
        print(f"  ❌ {sid}: No FOD template at {fod_template}")
        continue

    print(f"  {sid}:")
    ok = True
    for hemi, src in [('left', src_left), ('right', src_right)]:
        dst = cst_dir / f"roi_CST_{hemi}_fod.mif"
        cmd = (f"mrgrid {src} regrid {dst} "
               f"-template {fod_template} -interp nearest -force -quiet")
        print(f"    {hemi}: {src.parent.name}/{src.name} → {dst.name} ... ", end="", flush=True)
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
        if result.returncode == 0 and dst.exists():
            print("✅")
        else:
            print(f"❌ {result.stderr[:200]}")
            ok = False
    if ok:
        regrid_ok += 1

print(f"\n  Regridded {regrid_ok}/{len(all_to_regrid)} subjects successfully")

# %% [Step 5e-CST-QC] QC: Verify all AP/APTR subjects have CST ROIs in FOD space
print(f"{'=' * 70}")
print(f" CST ROI FINAL STATUS — ALL AP/APTR SUBJECTS")
print(f"{'=' * 70}")
all_ok = True
for s in subjects:
    if s['group'] == 'CTRL':
        continue
    sid = s['id']
    manual_roi_dir = s['base'] / 'manual_roi'
    cst_dir = s['wb_dir'] / "passthrough_cst"
    src_left  = manual_roi_dir / "roi_CST_left_flip.mif"
    src_right = manual_roi_dir / "roi_CST_right_flip.mif"
    fod_left  = cst_dir / "roi_CST_left_fod.mif"
    fod_right = cst_dir / "roi_CST_right_fod.mif"
    has_src  = src_left.exists() and src_right.exists()
    has_fod  = fod_left.exists() and fod_right.exists()
    icon = "✅" if has_fod else ("⚠️ " if has_src else "❌")
    if not has_fod:
        all_ok = False
    print(f"  {icon} [{s['group']:>4}] {sid:<30} "
          f"src={'YES' if has_src else 'NO':>3}  "
          f"fod_L={'YES' if fod_left.exists() else 'NO':>3}  "
          f"fod_R={'YES' if fod_right.exists() else 'NO':>3}")
if all_ok:
    print("\n  ✅ All AP/APTR subjects have CST ROIs — ready for Step 5f")
else:
    print("\n  ⚠️  Some subjects still missing CST ROIs — complete drawing before Step 5f")

# %% [Step 5e-CST-extract] Extract CST Passthrough Tracts (AP/APTR only) — PARALLEL
# ────────────────────────────────────────────────────────────────────────
# For each AP/APTR subject, extract streamlines that pass through BOTH:
#   1. The CST ROI  (passthrough_cst/roi_CST_{hemi}_fod.mif)
#   2. The WHOLE graft ROI (manual_roi/fod_space/roi_Graft_{hemi}_flip_fod.mif)
#
# Uses tckedit with dual -include from the whole-brain tractogram.
# Also extracts using shrunk graft ROI if available.
# Output → passthrough_cst/graft_{hemi}_to_CST_{hemi}.tck
# ────────────────────────────────────────────────────────────────────────
print(f"{'=' * 70}")
print(f" EXTRACTING CST PASSTHROUGH TRACTS (AP/APTR only)")
print(f"{'=' * 70}")

HEMIS = ['left', 'right']

# First, delete old segmented CST tracts (from previous dorsal-segment approach)
for s in subjects:
    if s['group'] == 'CTRL':
        continue
    cst_dir = s['wb_dir'] / "passthrough_cst"
    if not cst_dir.exists():
        continue
    old_tcks = [t for t in cst_dir.glob('seg9_graft_*.tck')]
    if old_tcks:
        for t in old_tcks:
            t.unlink()
        print(f"  🗑  {s['id']}: deleted {len(old_tcks)} old segmented CST tracts")

cst_tasks = []
for s in subjects:
    if s['group'] == 'CTRL':
        continue
    sid = s['id']
    wb_tck = s['wb_tck']
    cst_dir = s['wb_dir'] / "passthrough_cst"
    graft_fod_dir = s['base'] / 'manual_roi' / 'fod_space'
    graft_shrink_dir = graft_fod_dir / 'shrink'

    if not wb_tck.exists():
        print(f"  ❌ {sid}: no whole-brain tractogram, skipping")
        continue

    cmds = []
    for hemi in HEMIS:
        cst_roi = cst_dir / f"roi_CST_{hemi}_fod.mif"
        if not cst_roi.exists():
            print(f"  ⚠️  {sid}: no CST ROI for {hemi}, skipping")
            continue
        # --- Normal whole graft ---
        graft_roi = graft_fod_dir / f"roi_Graft_{hemi}_flip_fod.mif"
        out_tck = cst_dir / f"graft_{hemi}_to_CST_{hemi}.tck"
        if graft_roi.exists() and not out_tck.exists():
            cmds.append(
                f"{MODULES} && tckedit {wb_tck} {out_tck} "
                f"-include {graft_roi} -include {cst_roi} -force -quiet"
            )
        elif not graft_roi.exists():
            print(f"  ⚠️  {sid}: no whole graft ROI for {hemi}")
        # --- Shrunk whole graft ---
        graft_shrink_roi = graft_shrink_dir / f"roi_Graft_{hemi}_flip_fod_shrink.mif"
        out_shrink_tck = cst_dir / f"graft_{hemi}_to_CST_{hemi}_shrink.tck"
        if graft_shrink_roi.exists() and not out_shrink_tck.exists():
            cmds.append(
                f"{MODULES} && tckedit {wb_tck} {out_shrink_tck} "
                f"-include {graft_shrink_roi} -include {cst_roi} -force -quiet"
            )

    if cmds:
        compound = ' && '.join(cmds)
        cst_tasks.append((sid, compound, None))
    else:
        print(f"  ⚪ {sid}: all CST tracts already exist or no graft ROIs")

if cst_tasks:
    cst_results = run_parallel(cst_tasks, label="CST_passthrough")
else:
    print("  Nothing to extract.")

# Show results
print(f"\n  CST tract counts per subject:")
for s in subjects:
    if s['group'] == 'CTRL':
        continue
    cst_dir = s['wb_dir'] / "passthrough_cst"
    tcks = list(cst_dir.glob('*.tck'))
    n_normal = len([t for t in tcks if '_shrink' not in t.name and '_cleaned' not in t.name])
    n_shrink = len([t for t in tcks if '_shrink' in t.name and '_cleaned' not in t.name])
    print(f"  {s['id']:<30}: {n_normal} normal, {n_shrink} shrink")

# %% [Step 5e-CST-tckgen] Targeted CST tractography: seed CST → include Graft (AP/APTR) — PARALLEL
# ────────────────────────────────────────────────────────────────────────
# Runs tckgen seeded from the CST ROI, with the whole graft as an -include ROI.
# Uses the same iFOD2 parameters as the whole-brain tractogram (step=0.02,
# angle=60, cutoff=0.1, minlength=2, maxlength=15) but with -select 10000
# to get a reasonable number of targeted streamlines.
#
# Output → passthrough_cst/tckgen_CST_{hemi}_to_graft_{hemi}.tck
# ────────────────────────────────────────────────────────────────────────
print(f"{'=' * 70}")
print(f" TARGETED CST TRACTOGRAPHY: seed CST → include Graft (AP/APTR)")
print(f"{'=' * 70}")

HEMIS = ['left', 'right']
CST_TCKGEN_SELECT = 10000  # number of streamlines to select

tckgen_tasks = []
for s in subjects:
    if s['group'] == 'CTRL':
        continue
    sid = s['id']
    fod = s['dwi_dir'] / "wm_fod_norm_group.mif"
    mask = s['dwi_dir'] / "mask_anat.mif"
    cst_dir = s['wb_dir'] / "passthrough_cst"
    graft_fod_dir = s['base'] / 'manual_roi' / 'fod_space'

    if not fod.exists():
        print(f"  ❌ {sid}: no FOD image, skipping")
        continue
    if not mask.exists():
        print(f"  ❌ {sid}: no mask, skipping")
        continue

    cmds = []
    for hemi in HEMIS:
        cst_roi = cst_dir / f"roi_CST_{hemi}_fod.mif"
        graft_roi = graft_fod_dir / f"roi_Graft_{hemi}_flip_fod.mif"
        out_tck = cst_dir / f"tckgen_CST_{hemi}_to_graft_{hemi}.tck"

        if not cst_roi.exists():
            print(f"  ⚠️  {sid}: no CST ROI for {hemi}, skipping")
            continue
        if not graft_roi.exists():
            print(f"  ⚠️  {sid}: no Graft ROI for {hemi}, skipping")
            continue
        if out_tck.exists():
            continue

        cmds.append(
            f"{MODULES} && tckgen {fod} {out_tck} "
            f"-algorithm iFOD2 "
            f"-seed_image {cst_roi} "
            f"-include {graft_roi} "
            f"-mask {mask} "
            f"-select {CST_TCKGEN_SELECT} "
            f"-step 0.02 -angle 60 -cutoff 0.05 "
            f"-minlength 2 -maxlength 30 "
            f"-nthreads 4 -force -quiet"
        )

    if cmds:
        compound = ' && '.join(cmds)
        tckgen_tasks.append((sid, compound, None))
    else:
        print(f"  ⚪ {sid}: all tckgen CST tracts already exist or missing ROIs")

if tckgen_tasks:
    tckgen_results = run_parallel(tckgen_tasks, label="CST_tckgen")
else:
    print("  Nothing to generate.")

# Show results
print(f"\n  tckgen CST tract counts per subject:")
for s in subjects:
    if s['group'] == 'CTRL':
        continue
    cst_dir = s['wb_dir'] / "passthrough_cst"
    for hemi in HEMIS:
        out_tck = cst_dir / f"tckgen_CST_{hemi}_to_graft_{hemi}.tck"
        if out_tck.exists():
            import subprocess as _sp
            r = _sp.run(f"tckinfo {out_tck} 2>/dev/null | awk '/^ *count:/{{print $2}}'",
                        shell=True, capture_output=True, text=True)
            cnt = r.stdout.strip() or '?'
            print(f"  {s['id']:<30} {hemi:<6}: {cnt} streamlines")
        else:
            print(f"  {s['id']:<30} {hemi:<6}: NOT FOUND")

# %% [TEMP] CST tckgen — exclude Striatum/Cerebellum/Paleocortex/Thalamus (Apallial17NmF only)
# ────────────────────────────────────────────────────────────────────────
# Takes the tckgen CST→Graft tracts and removes any streamline passing
# through Striatum, Cerebellum, Paleocortex or Thalamus (both hemispheres).
# Output saved to /tmp/cst_viz/ for quick visualization — not stored officially.
# ────────────────────────────────────────────────────────────────────────
import subprocess

sid = "sub-Apallial17NmF"
ses = "ses-240816091633"
base = DERIV_BASE / sid / ses
roi_dir = base / "manual_roi" / "fod_space"
cst_dir = base / "dwi" / "wholebrain" / "passthrough_cst"

# Exclusion ROIs (all 8: 4 regions × 2 hemispheres)
exclude_rois = []
for region in ["Striatum", "cerebellum", "Palecortex", "Thalamus"]:
    for hemi in ["left", "right"]:
        roi = roi_dir / f"roi_{region}_{hemi}_flip_fod.mif"
        if roi.exists():
            exclude_rois.append(str(roi))
            print(f"  ✅ Exclude: {roi.name}")
        else:
            print(f"  ⚠️  Missing: {roi.name}")

exclude_flags = " ".join(f"-exclude {r}" for r in exclude_rois)
out_dir = Path("/tmp/cst_viz")
out_dir.mkdir(exist_ok=True)

for hemi in ["left", "right"]:
    src = cst_dir / f"tckgen_CST_{hemi}_to_graft_{hemi}.tck"
    dst = out_dir / f"Apallial17_CST_{hemi}_noTargets.tck"
    if not src.exists():
        print(f"  ❌ Source not found: {src.name}")
        continue
    cmd = f"{MODULES} && tckedit {src} {dst} {exclude_flags} -force -quiet"
    print(f"\n  Processing {hemi}...")
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
    if r.returncode == 0:
        # Count remaining streamlines
        info = subprocess.run(f"tckinfo {dst} 2>/dev/null | awk '/^ *count:/{{print $2}}'",
                              shell=True, capture_output=True, text=True)
        cnt = info.stdout.strip() or '?'
        print(f"  ✅ {dst.name}: {cnt} streamlines remaining")
    else:
        print(f"  ❌ Error: {r.stderr[:200]}")

print(f"\n  Output in: {out_dir}/")

# %% [Step 5e-CST-clean] Clean CST Passthrough Tracts with MATLAB AFQ — PARALLEL
# Same MATLAB AFQ cleaning as Step 5g, but for CST tracts in passthrough_cst/
print(f"{'=' * 70}")
print(f" CLEANING CST PASSTHROUGH TRACTS WITH MATLAB AFQ")
print(f"{'=' * 70}")

cst_clean_tasks = []
for s in subjects:
    if s['group'] == 'CTRL':
        continue
    sid = s['id']
    cst_dir = s['wb_dir'] / "passthrough_cst"
    raw_tcks = [t for t in cst_dir.glob('*.tck')
                if '_cleaned' not in t.name and '_SF' not in t.name]
    # Check if any need cleaning (no matching _cleaned file)
    needs_clean = [t for t in raw_tcks
                   if not (t.parent / t.name.replace('.tck', '_cleaned.tck')).exists()]
    if not needs_clean:
        print(f"  ⚪ {sid}: all CST tracts already cleaned")
        continue

    # Use clean_tracts.py with --passthrough-cst flag if it exists,
    # otherwise build MATLAB commands directly
    base_cmd = (f"{MODULES} && {MODULE_LOAD_MATLAB} && cd {REPO_DIR} && python "
                f"functions/clean_tracts.py --subject {sid} --passthrough-cst "
                f"--derivatives-dir {DERIV_BASE}")
    cst_clean_tasks.append((sid, base_cmd, None))

if cst_clean_tasks:
    cst_clean_results = run_parallel(cst_clean_tasks, label="CST_cleaning")
else:
    print("  Nothing to clean (all already done or no raw CST tracts).")

# Show cleaned counts
for s in subjects:
    if s['group'] == 'CTRL':
        continue
    cst_dir = s['wb_dir'] / "passthrough_cst"
    n_cln = len(list(cst_dir.glob('*_cleaned.tck')))
    print(f"  {s['id']:<30}: {n_cln} cleaned CST tracts")

# %% [markdown] Step 5f: Passthrough Tract Extraction
# ---
# ## Step 5f: Passthrough Tract Extraction
#
# Extract streamlines that pass through specific ROI pairs from the whole-brain tractogram.
# Uses `tckedit` with `-include` flags.
#
# Parameters from `config1.yaml → passthrough_extraction`.
# For AP/APTR: also extracts using **shrunk** (eroded) graft segments.

# %% [Step 5f] Show Passthrough Extraction Parameters
print("Passthrough extraction parameters (from config1.yaml):")
print(f"  Targets CTRL:    {pt.get('targets', {}).get('CTRL', [])}")
print(f"  Targets AP/APTR: {pt.get('targets', {}).get('AP_APTR', [])}")
print(f"  Segments CTRL:   {pt.get('segments', {}).get('CTRL', [])}")
print(f"  Segments AP:     {pt.get('segments', {}).get('AP_APTR', [])}")
print(f"  Hemispheres:     {pt.get('hemispheres', [])}")
print(f"  Cerebellum incl thalamus: {pt.get('cerebellum_include_thalamus', False)}")
print(f"  Shrink enabled:  {sk.get('enabled', False)}")
print(f"  Erosion voxels:  {sk.get('erosion_voxels', 'N/A')}")
print(f"  SEG_MODE:        {SEG_MODE}")
print(f"  To change: edit config1.yaml → passthrough_extraction / shrink_analysis")

# %% [Step 5f] Extract Passthrough Tracts — PARALLEL (both hemispheres + shrink)
# IMPORTANT: --create-shrink-segments processes ALL subjects globally, so it must
# run ONCE sequentially before the parallel extraction to avoid race conditions.
print(f"{'=' * 70}")
print(f" EXTRACTING PASSTHROUGH TRACTS (PARALLEL) — {len(subjects)} SUBJECTS")
print(f"{'=' * 70}")

# --- Step 1: Create shrunk segments ONCE (sequential, affects all AP/APTR) ---
has_ap_aptr = any(s['group'] in ['AP', 'APTR'] for s in subjects)
if has_ap_aptr:
    print("  Creating shrunk graft segments (sequential, all AP/APTR subjects)...")
    first_ap = next(s for s in subjects if s['group'] in ['AP', 'APTR'])
    shrink_cmd = (f"{MODULES} && cd {REPO_DIR} && python "
                  f"functions/extract_passthrough_tracts.py --subject {first_ap['id']} "
                  f"--create-shrink-segments --derivatives-dir {DERIV_BASE}")
    shrink_result = subprocess.run(shrink_cmd, shell=True, capture_output=True, text=True)
    if shrink_result.returncode == 0:
        print("  ✅ Shrunk segments created successfully")
    else:
        print(f"  ❌ Shrunk segments FAILED (rc={shrink_result.returncode})")
        print(f"     {shrink_result.stderr[-300:]}")

# --- Step 2: Extract passthrough tracts in parallel (left + right + shrink) ---
pt_tasks = []
for s in subjects:
    cmds = []
    base_cmd = (f"{MODULES} && cd {REPO_DIR} && python "
                f"functions/extract_passthrough_tracts.py --subject {s['id']}")
    cmds.append(f"{base_cmd} --hemi left --derivatives-dir {DERIV_BASE}")
    cmds.append(f"{base_cmd} --hemi right --derivatives-dir {DERIV_BASE}")
    if s['group'] in ['AP', 'APTR']:
        cmds.append(f"{base_cmd} --hemi left --shrink --derivatives-dir {DERIV_BASE}")
        cmds.append(f"{base_cmd} --hemi right --shrink --derivatives-dir {DERIV_BASE}")
    compound = ' && '.join(cmds)
    pt_tasks.append((s['id'], compound, None))

pt_results = run_parallel(pt_tasks, label="Passthrough extraction")

# %% [Step 5f-QC] QC: List tract counts per subject
print(f"{'=' * 70}")
print(f" PASSTHROUGH TRACT COUNTS")
print(f"{'=' * 70}")
for s in subjects:
    pt_d = s['pt_dir']
    if pt_d.exists():
        tracts = sorted(glob.glob(str(pt_d / '*.tck')))
        n_raw = len([t for t in tracts if '_cleaned' not in t])
        print(f"  {s['id']:<30} ({s['group']:<4}): {n_raw} raw tracts")
    else:
        print(f"  {s['id']:<30} ({s['group']:<4}): NO passthrough dir")

# %% [Step 5f-QC2] QC: View sample tract for first subject
s = subjects[0]
pt_d = s['pt_dir']
if pt_d.exists():
    sample = list(sorted(pt_d.glob('*thalamus_left*.tck')))
    sample = sample[0] if sample else None
    if sample:
        if s['group'] == 'CTRL':
            roi_overlay = s['base'] / 'atlas' / 'rois' / 'fod_space' / 'roi_thalamus_left.nii.gz'
        else:
            roi_overlay = s['base'] / 'manual_roi' / 'fod_space' / 'roi_Thalamus_left_flip_fod.mif'
        launch_gui(
            f"mrview {s['t2w']} "
            f"-overlay.load {roi_overlay} -overlay.opacity 0.4 "
            f"-tractography.load {sample}"
        )
        print(f"Showing: {sample.name} for {s['id']}")
    else:
        print("No thalamus tract found for visualization")

# %% [markdown] Step 5g: MATLAB AFQ Cleaning
# ---
# ## Step 5g: MATLAB AFQ Cleaning
#
# Remove outlier streamlines using MATLAB `AFQ_removeFiberOutliers`.
# All parameters from `config1.yaml → matlab_cleaning`.

# %% [Step 5g] Show MATLAB AFQ Cleaning Parameters
print("MATLAB AFQ Cleaning Parameters (from config1.yaml → matlab_cleaning):")
print(f"  max_dist:        {mat['max_dist']}  ← Gaussian distance from core (SD)")
print(f"  max_len:         {mat['max_len']}  ← Length deviation (SD)")
print(f"  num_nodes:       {mat['num_nodes']}  ← Resampling nodes per fiber")
print(f"  min_streamlines: {mat['min_streamlines']}  ← Skip cleaning if fewer")
print(f"  max_iter:        {mat.get('max_iter', 'N/A')}  ← Max cleaning iterations")
print(f"  core_method:     {mat.get('core_method', 'N/A')}")
print(f"  matlab_bin:      {mat.get('matlab_bin', 'N/A')}")
print(f"  afq_path:        {mat.get('afq_path', 'N/A')}")
print(f"  To change: edit config1.yaml → matlab_cleaning")

# %% [Step 5g] Run MATLAB AFQ Cleaning — PARALLEL
# NOTE: MATLAB uses ~2 cores per instance => 12 subjects × 2 = 24 cores (fine for lmc02)
print(f"{'=' * 70}")
print(f" CLEANING PASSTHROUGH TRACTS WITH MATLAB AFQ (PARALLEL)")
print(f"{'=' * 70}")

clean_tasks = []
for s in subjects:
    base_cmd = (f"{MODULES} && cd {REPO_DIR} && python "
                f"functions/clean_tracts.py --subject {s['id']} --passthrough "
                f"--derivatives-dir {DERIV_BASE}")
    if s['group'] in ['AP', 'APTR']:
        # Clean both normal and shrunk tracts
        compound = f"{base_cmd} && {base_cmd} --shrink"
        clean_tasks.append((s['id'], compound, 7200))
    else:
        clean_tasks.append((s['id'], base_cmd, 3600))

clean_results = run_parallel(clean_tasks, label="MATLAB cleaning")

# %% [Step 5g-QC] QC: Compare raw vs cleaned for first subject
s = subjects[0]
cleaned_files = sorted(glob.glob(str(s['pt_dir'] / '*_cleaned.tck'))) if s['pt_dir'].exists() else []
if cleaned_files:
    cleaned = Path(cleaned_files[0])
    raw = cleaned.parent / cleaned.name.replace('_cleaned.tck', '.tck')
    print(f"Comparing: {raw.name}")
    run(f"{MODULES} && tckinfo {raw} -count")
    print(f"Cleaned: {cleaned.name}")
    run(f"{MODULES} && tckinfo {cleaned} -count")
    launch_gui(f"mrview {s['t2w']} -tractography.load {cleaned}")
    print("Cleaned tract should be a tight, clean bundle.")
else:
    print("No cleaned tracts found yet.")

# %% [Step 5fg-CTRL-rerun] Re-run passthrough + cleaning for CTRL only (after ROI L/R fix)
# ────────────────────────────────────────────────────────────────────────
# After fixing isocortex & olfactory L/R labels for CTRL subjects, we must:
#   1. Delete old passthrough tracts (raw + cleaned + SF)
#   2. Re-extract passthrough tracts  (Step 5f logic, CTRL only)
#   3. Re-run MATLAB AFQ cleaning     (Step 5g logic, CTRL only)
#
# AP and APTR are NOT affected — their ROIs were verified correct.
# ────────────────────────────────────────────────────────────────────────
import shutil

ctrl_subjects = [s for s in subjects if s['group'] == 'CTRL']
print(f"{'=' * 70}")
print(f" CTRL-ONLY RERUN: PASSTHROUGH + CLEANING ({len(ctrl_subjects)} subjects)")
print(f"{'=' * 70}")

# ── Part 1: Delete old CTRL passthrough tracts ──
print(f"\n{'─' * 50}")
print(f" Part 1: Deleting old CTRL passthrough tracts")
print(f"{'─' * 50}")
for s in ctrl_subjects:
    pt_d = s['pt_dir']  # .../wholebrain/passthrough_tracts/
    if pt_d.exists():
        old_tcks = list(pt_d.glob('*.tck'))
        if old_tcks:
            for f in old_tcks:
                f.unlink()
            print(f"  🗑  {s['id']}: deleted {len(old_tcks)} old .tck files")
        else:
            print(f"  ⚪ {s['id']}: no .tck files to delete")
    else:
        print(f"  ⚪ {s['id']}: passthrough_tracts/ does not exist")

# ── Part 2: Re-extract passthrough tracts (CTRL only, left + right) ──
print(f"\n{'─' * 50}")
print(f" Part 2: Extracting passthrough tracts (CTRL × left + right)")
print(f"{'─' * 50}")

pt_ctrl_tasks = []
for s in ctrl_subjects:
    base_cmd = (f"{MODULES} && cd {REPO_DIR} && python "
                f"functions/extract_passthrough_tracts.py --subject {s['id']}")
    cmds = [
        f"{base_cmd} --hemi left --derivatives-dir {DERIV_BASE}",
        f"{base_cmd} --hemi right --derivatives-dir {DERIV_BASE}",
    ]
    compound = ' && '.join(cmds)
    pt_ctrl_tasks.append((s['id'], compound, None))

pt_ctrl_results = run_parallel(pt_ctrl_tasks, label="CTRL_passthrough")

# Show tract counts after extraction
for s in ctrl_subjects:
    pt_d = s['pt_dir']
    if pt_d.exists():
        n = len(list(pt_d.glob('*.tck')))
        print(f"  {s['id']}: {n} tracts extracted")

# ── Part 3: Re-run MATLAB AFQ cleaning (CTRL only) ──
print(f"\n{'─' * 50}")
print(f" Part 3: MATLAB AFQ cleaning (CTRL only)")
print(f"{'─' * 50}")

clean_ctrl_tasks = []
for s in ctrl_subjects:
    base_cmd = (f"{MODULES} && {MODULE_LOAD_MATLAB} && cd {REPO_DIR} && python "
                f"functions/clean_tracts.py --subject {s['id']} --passthrough "
                f"--derivatives-dir {DERIV_BASE}")
    clean_ctrl_tasks.append((s['id'], base_cmd, None))

clean_ctrl_results = run_parallel(clean_ctrl_tasks, label="CTRL_cleaning")

# Show cleaned counts
for s in ctrl_subjects:
    pt_d = s['pt_dir']
    if pt_d.exists():
        n_raw = len([f for f in pt_d.glob('*.tck') if '_cleaned' not in f.name])
        n_cln = len(list(pt_d.glob('*_cleaned.tck')))
        print(f"  {s['id']}: {n_raw} raw, {n_cln} cleaned")

print(f"\n{'=' * 70}")
print(f" CTRL RERUN COMPLETE")
print(f"{'=' * 70}")

# %% [markdown] Step 5h: Generate Metrics Report
# ---
# ## Step 5h: Generate Tract Metrics Report (v2)
#
# Compute tract metrics and generate CSV + HTML report:
# - **Streamline count** (from cleaned passthrough tracts)
# - **SIFT2 weighted** (sum of wholebrain SIFT2 weights for matched raw streamlines)
# - **Density** (streamline count / target ROI volume in mm³)
#
# Source ROIs: CTRL → isocortex {ant, cen, pos}; AP/APTR → seg9 graft dorsal {ant, cen, pos} (shrunk)
#
# Target ROIs: L/R × {Striatum, Thalamus, Paleocortex, Cerebellum}
#
# Output: CSV (main + 3 companion) + HTML (tables + bar charts with error bars + individual points)

# %% [Step 5h] Generate Tract Metrics Report v2 (CSV + HTML)
import datetime as _dt
report_version = f"v{_dt.datetime.now().strftime('%Y%m%d')}_newmask"

print(f"{'=' * 70}")
print(f" GENERATING TRACT METRICS REPORT v2")
print(f" Version: {report_version}")
print(f"{'=' * 70}")

report_cmd = (
    f"{MODULES} && cd {REPO_DIR} && "
    f"python functions/generate_tract_metrics_report_v2.py "
    f"--version {report_version} "
    f"--output-dir {REPO_DIR / 'reports'} "
    f"--derivatives-dir {DERIV_BASE}"
)
print(f"Command:\n  {report_cmd[:300]}...")
result = run(report_cmd, timeout=7200)

reports_dir = REPO_DIR / 'reports'
check_file(reports_dir / f'tract_metrics_{report_version}.csv', 'Main CSV')
check_file(reports_dir / f'tract_metrics_count_{report_version}.csv', 'Count CSV')
check_file(reports_dir / f'tract_metrics_sift2_{report_version}.csv', 'SIFT2 CSV')
check_file(reports_dir / f'tract_metrics_density_{report_version}.csv', 'Density CSV')
check_file(reports_dir / f'tract_metrics_report_{report_version}.html', 'HTML report')

# %% [Step 5h] Preview CSV Report
import pandas as pd
csv_path = REPO_DIR / "reports" / f"tract_metrics_{report_version}.csv"
if csv_path.exists():
    df = pd.read_csv(csv_path)
    print(f"Report: {len(df)} rows\n")
    print("Columns:", list(df.columns))
    print(f"\nFirst 20 rows:")
    print(df.head(20).to_string(index=False))
else:
    print("Report CSV not found.")

# %% [markdown] Step 5i: Fixel-Based Analysis (FBA)
# ---
# ## Step 5i: Fixel-Based Analysis (FBA)
#
# Decompose FODs into discrete fixels and compute microstructural metrics.
#
# **FBA pipeline:**
# 1. **Generate fixels** — `fod2fixel` on `wm_fod_norm_group.mif` → `fixels_group/`
# 2. **Voxel AFD map** — `fixel2voxel` collapses fixel AFD to voxel space → `afd_voxel_sum.mif`
# 3. **ROI fixel metrics** — directional FD within each ROI (AP/DV/LR decomposition)
# 4. **Tract AFD sampling** — `tcksample` samples AFD along cleaned passthrough tracts
#
# **Key metrics:**
# - `fd_total` — Total Fiber Density in ROI
# - `fd_ap`, `fd_dv`, `fd_lm` — Directional FD (Anterior-Posterior, Dorsal-Ventral, Lateral-Medial)
# - `mean_afd` — Mean AFD along a tract path
#
# **Axis convention (after flip1):** axis 0 = AP, axis 1 = DV, axis 2 = LR

# %% [Step 5i-a] Generate Fixels (fod2fixel) — PARALLEL
print(f"{'=' * 70}")
print(f" GENERATING FIXELS (fod2fixel) — PARALLEL — {len(subjects)} SUBJECTS")
print(f"{'=' * 70}")

fixel_tasks = []
for s in subjects:
    DWI = s['dwi_dir']
    FIXEL_DIR = s['fixel_dir']
    fod_file = DWI / 'wm_fod_norm_group.mif'
    mask_file = DWI / 'mask_anat.mif'
    if not fod_file.exists():
        print(f"  ❌ {s['id']}: no group FOD, skipping")
        continue
    # fod2fixel requires an empty output dir — remove stale fixels first
    fixel_tasks.append((
        s['id'],
        f"{MODULES} && rm -rf {FIXEL_DIR} && "
        f"fod2fixel {fod_file} {FIXEL_DIR} "
        f"-mask {mask_file} -afd afd.mif -peak_amp peak.mif -force -quiet",
        600
    ))

fixel_results = run_parallel(fixel_tasks, label="Fixel generation")

# %% [Step 5i-b] Create Voxel AFD Maps (fixel2voxel) — PARALLEL
print(f"{'=' * 70}")
print(f" CREATING VOXEL AFD MAPS — PARALLEL — {len(subjects)} SUBJECTS")
print(f"{'=' * 70}")

afd_tasks = []
for s in subjects:
    FIXEL_DIR = s['fixel_dir']
    AFD_VOXEL = s['afd_voxel']
    if (FIXEL_DIR / 'afd.mif').exists():
        afd_tasks.append((
            s['id'],
            f"{MODULES} && fixel2voxel {FIXEL_DIR}/afd.mif sum {AFD_VOXEL} -force",
            120
        ))
    else:
        print(f"  ❌ {s['id']}: no fixels, skipping")

afd_voxel_results = run_parallel(afd_tasks, label="Voxel AFD maps")

# %% [Step 5i-c] ROI-Based Fixel Metrics (Directional FD) — PARALLEL
print(f"{'=' * 70}")
print(f" COMPUTING ROI FIXEL METRICS — PARALLEL — {len(subjects)} SUBJECTS")
print(f"{'=' * 70}")

roi_fixel_tasks = []
for s in subjects:
    if s['group'] == 'CTRL':
        roi_names = ['isocortex_left', 'isocortex_right',
                     'thalamus_left', 'thalamus_right',
                     'striatum_left', 'striatum_right',
                     'cerebellum_left', 'cerebellum_right']
    else:
        roi_names = ['graft_left', 'graft_right',
                     'thalamus_left', 'thalamus_right',
                     'striatum_left', 'striatum_right',
                     'palecortex_left', 'palecortex_right',
                     'cerebellum_left', 'cerebellum_right']
    fixel_csv = REPO_DIR / 'reports' / f'fixel_metrics_{s["id"]}.csv'
    # Chain all ROIs for this subject into one compound command
    cmds = [f"{MODULES} && cd {REPO_DIR} && python functions/analyze_roi_fixels.py "
            f"--subject {s['id']} --roi {roi} --output {fixel_csv} --force --derivatives-dir {DERIV_BASE}"
            for roi in roi_names]
    roi_fixel_tasks.append((s['id'], ' && '.join(cmds), 1200))

roi_fixel_results = run_parallel(roi_fixel_tasks, label="ROI fixel metrics")

# %% [Step 5i-d] Tract AFD Sampling (tcksample) — PARALLEL
import numpy as np
print(f"{'=' * 70}")
print(f" SAMPLING AFD ALONG CLEANED TRACTS — PARALLEL — {len(subjects)} SUBJECTS")
print(f"{'=' * 70}")

tcksample_tasks = []
for s in subjects:
    AFD_VOXEL = s['afd_voxel']
    if not AFD_VOXEL.exists():
        print(f"  ❌ {s['id']}: no afd_voxel_sum.mif, skipping")
        continue
    # Collect all cleaned tracts for this subject
    cleaned = sorted(glob.glob(str(s['pt_dir'] / '*_cleaned.tck'))) if s['pt_dir'].exists() else []
    if s['group'] in ['AP', 'APTR']:
        cleaned += sorted(glob.glob(str(s['pt_shrink_dir'] / '*_cleaned.tck'))) if s['pt_shrink_dir'].exists() else []
    if not cleaned:
        print(f"  ⚠️  {s['id']}: no cleaned tracts")
        continue
    # Build compound tcksample command for all tracts of this subject
    cmds = []
    for tp in cleaned:
        tn = Path(tp).stem.replace('_cleaned', '')
        afd_out = Path(tp).parent / f'{tn}_afd.txt'
        cmds.append(f"{MODULES} && tcksample {tp} {AFD_VOXEL} {afd_out} -stat_tck mean -force")
    tcksample_tasks.append((s['id'], ' && '.join(cmds), 600))

tcksample_results = run_parallel(tcksample_tasks, label="Tract AFD sampling")

# %% [Step 5i-e] Slicewise Fixel Analysis + 3-Group Comparison
print(f"{'=' * 70}")
print(f" SLICEWISE FIXEL ANALYSIS — 3-GROUP COMPARISON")
print(f"{'=' * 70}")

# Run the 3-group comparison (compares CTRL vs AP vs APTR per slice)
for hemi in ['left', 'right']:
    print(f"\n  --- {hemi.upper()} hemisphere ---")
    run(f"{MODULES} && cd {REPO_DIR} && python functions/analyze_roi_fixels.py "
        f"--compare --force --derivatives-dir {DERIV_BASE}",
        timeout=600)

print("\n--- Slicewise fixel reports ---")
for hemi in ['left', 'right']:
    for shrink_tag in ['', '_shrink']:
        fig_path = REPO_DIR / 'reports' / f'fixel_3group_comparison_{hemi}{shrink_tag}.png'
        if fig_path.exists():
            check_file(fig_path)

# %% [Step 5i-km] 3-Group Fixel K-means Direction Comparison
import sys, yaml, importlib
import numpy as np
from pathlib import Path
sys.path.insert(0, str(REPO_DIR))

# Ensure derivatives override is set for direct Python imports
os.environ['MBRAIN_DERIVATIVES_DIR'] = str(DERIV_BASE)

import functions.fixel_roi_analysis as _fri
importlib.reload(_fri)
from functions.fixel_roi_analysis import analyze_3group_kmeans

print('Step 5i-km: 3-Group Fixel K-means Direction Comparison\n')

with open(REPO_DIR / 'config1.yaml') as f:
    cfg_km = yaml.safe_load(f)

# Override derivatives paths in the loaded config
cfg_km.setdefault('paths', {})['derivatives_flip1'] = str(DERIV_BASE)
cfg_km['paths']['derivatives_dir'] = str(DERIV_BASE)

report_dir = REPO_DIR / 'reports'
combos = [(h, True) for h in ['left', 'right']]

# Pass 1: compute data to find unified y_max
print('Pass 1: computing data for all combinations...')
all_results_km = {}
global_max = 0.0
for hemi, use_shrink in combos:
    tag = f'{hemi}_shrink'
    r = analyze_3group_kmeans(
        cfg_km, hemisphere=hemi, use_shrink=True,
        output_dir=report_dir, n_clusters=2, conda_env='mbrain',
        force=False, show=False, verbose=False, plot=False)
    all_results_km[tag] = (hemi, True, r)
    if r:
        for group in ['CTRL', 'AP', 'APTR']:
            for subj, data in r['all_group_data'].get(group, {}).items():
                for s_idx, res in data.items():
                    for sd in res.get('scaled_directions', []):
                        val = np.max(np.abs(sd))
                        if val > global_max:
                            global_max = val

y_max_unified = max(0.05, np.ceil(global_max * 20) / 20)
print(f'  Unified y_max = {y_max_unified:.4f} (global max = {global_max:.4f})\n')

# Pass 2: generate plots with forced y_max
for tag, (hemi, use_shrink, r) in all_results_km.items():
    print(f'\n{"="*60}')
    print(f'  {tag.upper()}')
    print(f'{"="*60}')
    if not r:
        print(f'  ⚠ No data for {tag}')
        continue
    from functions.fixel_roi_analysis import (
        plot_3group_kmeans_comparison, plot_3group_component_lines
    )
    suffix = '_shrink' if use_shrink else ''
    tns = r['target_n_slices']

    fig_path = report_dir / f'fixel_3group_comparison_{hemi}{suffix}.png'
    plot_3group_kmeans_comparison(
        r['all_group_data'], hemi, tns,
        use_shrink=use_shrink, save_path=fig_path, show=True,
        y_max=y_max_unified)

    lines_path = report_dir / f'fixel_3group_components_{hemi}{suffix}.png'
    plot_3group_component_lines(
        r['all_group_data'], hemi, tns,
        use_shrink=use_shrink, save_path=lines_path, show=True,
        y_max=y_max_unified)

    print(f'  ✅ Radar: {fig_path}')
    print(f'  ✅ Lines: {lines_path}')

# %% [Step 5i-fig] Generate Fixel 3×3 Composite Figure v08 (newmask)
# Generates the 3×3 composite figure:
#   Row 1: FOD screenshots (CTRL | AP | APTR)
#   Row 2: Fixel screenshots
#   Row 3: Per-mouse radar (4 mice × DV-dominant direction, first-thalamus slice)
# Uses mrtrix_newmask derivatives. Output: reports/fixel_4x3_comparison_v08.{png,svg}

print(f"{'=' * 70}")
print(f" GENERATING FIXEL 3×3 COMPOSITE FIGURE v08 [newmask]")
print(f"{'=' * 70}")

fig_cmd = (
    f"{MODULES} && cd {REPO_DIR} && "
    f"python tmp/make_4x3_figure_v08.py"
)
print(f"Command:\n  {fig_cmd}\n")
result = run(fig_cmd, timeout=3600)

check_file(REPO_DIR / 'reports' / 'fixel_4x3_comparison_v08.png', 'Fixel 4×3 PNG')
check_file(REPO_DIR / 'reports' / 'fixel_4x3_comparison_v08.svg', 'Fixel 4×3 SVG')

# %% [markdown] QC 5i: Fixel Analysis
# ### QC Checkpoint: Fixel Analysis
#
# **Verify:**
# - Fixels generated without errors
# - AFD voxel map looks reasonable (higher in WM, low in CSF)
# - ROI metrics show expected directional patterns
# - Tract AFD values consistent across subjects in same group

# %% [QC 5i] View AFD Voxel Map — check all subjects
print('QC: Checking AFD voxel maps for all subjects\n')
found_afd = 0
missing_afd = 0
for s in subjects:
    if s['afd_voxel'].exists():
        print(f"  ✅ {s['id']}: afd_voxel_sum.mif exists")
        found_afd += 1
    else:
        print(f"  ❌ {s['id']}: afd_voxel_sum.mif NOT FOUND")
        missing_afd += 1

print(f"\n  Summary: {found_afd} found, {missing_afd} missing")

# Open first available for visual check
for s in subjects:
    if s['afd_voxel'].exists():
        print(f"\nOpening {s['id']} for visual check...")
        launch_gui(f"mrview {s['t2w']} -overlay.load {s['afd_voxel']} -overlay.opacity 0.5")
        print("AFD voxel map: should be bright in WM regions, dark in CSF/outside brain.")
        break

# %% [markdown] Step 5j: Fixel Square ROI Analysis
# ---
# ## Step 5j: Fixel Square ROI Direction Analysis
#
# Analyze the **main fiber directions** within a user-drawn square ROI.
#
# **Pipeline:**
# 1. Run `fod2fixel` masked to the square ROI → local fixel directory
# 2. Extract fixel directions (unit vectors) and AFD values
# 3. Cluster fixels into 2 populations (K-means on direction vectors)
# 4. Report **Main Direction 1** and **Main Direction 2** with mean AFD
#
# **Prerequisites:**
# - `square_roi_fixel.mif` drawn in mrview ROI editor and saved in `dwi/`
# - `wm_fod_norm_group.mif` (from Step 3b / re-run in Step 5b)
#
# **Axis convention (flip1):** axis 0 = AP, axis 1 = DV, axis 2 = LR

# %% [Step 5j] Fixel Square ROI — Generate fixels for each subject
from functions.fixel_roi_analysis import generate_fixel_peaks

print(f"{'=' * 70}")
print(f" FIXEL SQUARE ROI ANALYSIS")
print(f"{'=' * 70}")

for i, s in enumerate(subjects):
    label = f"[{i+1}/{len(subjects)}] {s['id']} ({s['group']})"
    DWI = s['dwi_dir']
    SQUARE_ROI_T2W = DWI / 'square_roi_fixel.mif'
    FOD_PATH = DWI / 'wm_fod_norm_group.mif'
    
    if not SQUARE_ROI_T2W.exists():
        print(f"  ⏭️ {label}: no square_roi_fixel.mif, skipping")
        continue
    if not FOD_PATH.exists():
        print(f"  ❌ {label}: no wm_fod_norm_group.mif, skipping")
        continue
    
    peaks_nii = generate_fixel_peaks(
        roi_path=SQUARE_ROI_T2W,
        fod_path=FOD_PATH,
        output_dir=DWI,
        conda_env='mbrain',
        force=True,
    )
    print(f"  ✅ {label}: fixel extraction complete")

# %% [Step 5j-2] Fixel Square ROI — K-means clustering for each subject
import numpy as np
from functions.fixel_roi_analysis import extract_peak_vectors, cluster_fixel_directions, plot_fixel_directions

print(f"{'=' * 70}")
print(f" FIXEL SQUARE ROI — K-MEANS CLUSTERING")
print(f"{'=' * 70}")

for i, s in enumerate(subjects):
    label = f"[{i+1}/{len(subjects)}] {s['id']} ({s['group']})"
    DWI = s['dwi_dir']
    PEAKS_NII = DWI / 'square_roi_peaks.nii'
    
    if not PEAKS_NII.exists():
        print(f"  ⏭️ {label}: no peaks file, skipping")
        continue
    
    all_vectors, all_peak_indices, n_peaks = extract_peak_vectors(PEAKS_NII)
    if len(all_vectors) == 0:
        print(f"  ⚠️ {label}: no valid fixel vectors")
        continue
    
    print(f"\n  {label}: {len(all_vectors)} fixel vectors")
    cluster_result = cluster_fixel_directions(
        all_vectors, all_peak_indices, n_peaks,
        n_clusters=2, verbose=True,
    )
    
    if cluster_result is not None:
        fig_path = REPO_DIR / 'reports' / f'fixel_square_roi_{s["id"]}.png'
        fig = plot_fixel_directions(
            cluster_result,
            subject_id=s['id'],
            save_path=fig_path,
            show=True,
        )
        print(f"  ✅ {s['id']}: saved {fig_path.name}")

# %% [markdown] Step 6: QC — Verify New Masks in Context
# ---
# ## Step 6: Final QC — Verify New Masks + Downstream Outputs
#
# Quick visual checks after re-running the pipeline.

# %% [Step 6a] QC: View new mask on registered DWI for each subject
for i, s in enumerate(subjects):
    b0 = s['dwi_dir'] / "b0_final.mif"
    mask = s['dwi_dir'] / "mask_anat.mif"
    if b0.exists() and mask.exists():
        print(f"[{i+1}/{len(subjects)}] {s['id']} ({s['group']})")
        launch_gui(
            f"mrview {s['t2w']} "
            f"-overlay.load {mask} -overlay.opacity 0.3 -overlay.colour 0,1,0"
        )
        answer = input("  PASS (Enter) or FAIL (f): ").strip().lower()
        if answer == 'f':
            print(f"  ⚠️ {s['id']} FAILED final QC — may need re-editing")
    else:
        print(f"  ⏭️ {s['id']}: missing b0_final or mask_anat")

# %% [Step 6b] QC: Summary table — mask sizes
print(f"{'=' * 80}")
print(f" {'Subject':<30} {'Group':<6} {'Old mask MB':>12} {'New mask MB':>12} {'Ratio':>8}")
print(f"{'=' * 80}")

for s in subjects:
    old = s['dwi_dir'] / f"mask_old_{BACKUP_SUFFIX}.mif"
    new = s['dwi_dir'] / "mask.mif"
    
    old_sz = old.stat().st_size / (1024*1024) if old.exists() else 0
    new_sz = new.stat().st_size / (1024*1024) if new.exists() else 0
    ratio = new_sz / old_sz if old_sz > 0 else float('inf')
    
    flag = "" if 0.5 < ratio < 2.0 else " ⚠️"
    print(f" {s['id']:<30} {s['group']:<6} {old_sz:>10.3f}  {new_sz:>10.3f}  {ratio:>6.2f}x{flag}")

print(f"{'=' * 80}")

# %% [markdown] Step 7: Final Visualization & Summary
# ---
# ## Step 7: Final Visualization & Summary

# %% [Step 7a] Pipeline Summary Table — ALL subjects
print(f"\n{'=' * 80}")
print(f" PIPELINE COMPLETE: {len(subjects)} SUBJECTS")
print(f"{'=' * 80}\n")

for s in subjects:
    pt_d = s['pt_dir']
    raw_tracts = glob.glob(str(pt_d / '*.tck')) if pt_d.exists() else []
    raw_tracts = [t for t in raw_tracts if '_cleaned' not in t]
    cleaned_tracts = glob.glob(str(pt_d / '*_cleaned.tck')) if pt_d.exists() else []
    
    wb_ok = "✅" if s['wb_tck'].exists() else "❌"
    sift2_ok = "✅" if s['sift2_weights'].exists() else "❌"
    
    print(f"  {s['id']:<30} ({s['group']:<4})")
    print(f"    WB tractogram: {wb_ok}  |  SIFT2: {sift2_ok}")
    print(f"    Raw tracts: {len(raw_tracts):>3}  |  Cleaned: {len(cleaned_tracts):>3}")
    
    if s['group'] in ['AP', 'APTR']:
        shrink_d = s['pt_shrink_dir']
        shrink_raw = glob.glob(str(shrink_d / '*.tck')) if shrink_d.exists() else []
        shrink_raw = [t for t in shrink_raw if '_cleaned' not in t]
        shrink_cleaned = glob.glob(str(shrink_d / '*_cleaned.tck')) if shrink_d.exists() else []
        print(f"    Shrunk raw: {len(shrink_raw):>3}  |  Shrunk cleaned: {len(shrink_cleaned):>3}")

print(f"\nReport: reports/tract_metrics_{report_version}.csv")
print(f"HTML:   reports/tract_metrics_report_{report_version}.html")
print(f"\n{'=' * 80}")

# %% [Step 7b] Final Visualization — Cleaned Tracts for first subject per group
for group_name in ['CTRL', 'AP', 'APTR']:
    group_subjects = [s for s in subjects if s['group'] == group_name]
    if not group_subjects:
        continue
    s = group_subjects[0]
    pt_d = s['pt_dir']
    
    if s['group'] == 'CTRL':
        seg_roi = s['base'] / 'atlas' / 'rois' / 'fod_space' / 'roi_isocortex_left.nii.gz'
        target_roi = s['base'] / 'atlas' / 'rois' / 'fod_space' / 'roi_thalamus_left.nii.gz'
    else:
        seg_roi = s['base'] / 'manual_roi' / 'fod_space' / 'roi_Graft_left_flip_fod.mif'
        target_roi = s['base'] / 'manual_roi' / 'fod_space' / 'roi_Thalamus_left_flip_fod.mif'
    
    cleaned_left = sorted(glob.glob(str(pt_d / '*left*_cleaned.tck')))
    cmd_viz = (f"mrview {s['t2w']} "
               f"-overlay.load {seg_roi} -overlay.opacity 0.3 "
               f"-overlay.load {target_roi} -overlay.opacity 0.3")
    for t in cleaned_left[:6]:
        cmd_viz += f" -tractography.load {t}"
    
    print(f"Opening final visualization for {s['id']} ({group_name})...")
    launch_gui(cmd_viz)
    print(f"  Showing left-hemisphere cleaned tracts for {s['id']}")

# %% [markdown] Bar-plots: APTR SIFT2 Density
# ---
# ## Bar-plots: APTR SIFT2 Density (Left Hemisphere & Central Graft)
#
# Bar plots showing SIFT2 density (connectivity / seed_voxels × 1000) for all
# APTR mice. Individual data points overlaid on group means.

# %% [Bar-plots] APTR Density — Left hemisphere & Central summary
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Load metrics CSV (cleaned tracts)
csv_path = REPO_DIR / "reports" / f"tract_metrics_{report_version}.csv"
if not csv_path.exists():
    print(f"❌ Report CSV not found: {csv_path}")
    print("  Run Step 5h first to generate the report.")
else:
    df_plot = pd.read_csv(csv_path)
    
    # Filter APTR, left hemisphere only
    aptr_left = df_plot[(df_plot['group'] == 'APTR') & (df_plot['hemisphere'] == 'left')].copy()
    targets_order = ['olfactory_cortex', 'striatum', 'thalamus', 'cerebellum']
    target_labels = {'thalamus': 'Thalamus', 'striatum': 'Striatum',
                     'olfactory_cortex': 'Paleocortex', 'cerebellum': 'Cerebellum'}
    segments_order = ['anterior', 'central', 'posterior']
    seg_colors = {'anterior': '#E24A33', 'central': '#348ABD', 'posterior': '#988ED5'}
    
    plot_subjects = aptr_left['subject'].unique()
    n_subj = len(plot_subjects)
    marker_list = ['o', 's', 'D', '^']
    n_tgt = len(targets_order)
    
    def _single_seg_panel(ax, data_df, seg, color, title):
        """Draw single-segment bars for 4 targets, with individual points + subject legend."""
        seg_data = data_df[data_df['segment'] == seg]
        means, sems, all_vals = [], [], []
        for tgt in targets_order:
            vals = seg_data[seg_data['target'] == tgt]['sift2_density'].values
            means.append(np.mean(vals) if len(vals) > 0 else 0)
            sems.append(np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0)
            all_vals.append(vals)
        
        x_c = np.arange(n_tgt)
        ax.bar(x_c, means, 0.55, yerr=sems, capsize=4,
               color=color, alpha=0.55, edgecolor='black', linewidth=0.5, zorder=2)
        
        for tgt_idx in range(n_tgt):
            tgt_subjs = seg_data[seg_data['target'] == targets_order[tgt_idx]]['subject'].values
            for val_idx, val in enumerate(all_vals[tgt_idx]):
                subj_idx = list(plot_subjects).index(tgt_subjs[val_idx])
                short = tgt_subjs[val_idx].replace('sub-', '')
                jitter = (np.random.RandomState(subj_idx).rand() - 0.5) * 0.3
                label = short if tgt_idx == 0 else None
                ax.scatter(x_c[tgt_idx] + jitter, val,
                           marker=marker_list[subj_idx % len(marker_list)],
                           color=color, edgecolor='black', s=50, zorder=3,
                           linewidth=0.6, label=label)
        
        ax.set_xticks(x_c)
        ax.set_xticklabels([target_labels[t] for t in targets_order], fontsize=9,
                           rotation=30, ha='right')
        ax.set_title(title, fontsize=12)
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(title='Subject', fontsize=6, loc='upper right')
    
    import matplotlib.gridspec as gridspec
    from matplotlib.image import imread as mpl_imread
    
    fig = plt.figure(figsize=(20, 14))
    fig.suptitle('APTR — SIFT2 Density (Left Hemisphere, Cleaned Tracts)',
                 fontsize=14, fontweight='bold')
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3,
                           height_ratios=[1, 0.9])
    
    # Row 0: 3 segment bar subplots
    ax_ant = fig.add_subplot(gs[0, 0])
    _single_seg_panel(ax_ant, aptr_left, 'anterior', seg_colors['anterior'],
                      'Anterior Graft → Destinations')
    ax_ant.set_ylabel('SIFT2 Density\n(conn / seed_vox × 1000)', fontsize=10)
    
    ax_cen = fig.add_subplot(gs[0, 1])
    _single_seg_panel(ax_cen, aptr_left, 'central', seg_colors['central'],
                      'Central Graft → Destinations')
    ax_cen.set_ylabel('SIFT2 Density', fontsize=10)
    
    ax_pos = fig.add_subplot(gs[0, 2])
    _single_seg_panel(ax_pos, aptr_left, 'posterior', seg_colors['posterior'],
                      'Posterior Graft → Destinations')
    ax_pos.set_ylabel('SIFT2 Density', fontsize=10)
    
    # Row 1: mrview screenshots by destination
    destinations = ['paleocortex', 'striatum', 'thalamus']
    dest_labels = {'paleocortex': 'Paleocortex', 'striatum': 'Striatum', 'thalamus': 'Thalamus'}
    
    for col_idx, dest in enumerate(destinations):
        ax = fig.add_subplot(gs[1, col_idx])
        img_path = REPO_DIR / 'reports' / f'mrview_most_{dest}.png'
        title = f'{dest_labels[dest]}'
        if img_path.exists():
            img = mpl_imread(str(img_path))
            ax.imshow(img)
            ax.set_title(title, fontsize=10, fontweight='bold')
        else:
            ax.text(0.5, 0.5, f'Screenshot not yet saved:\n{img_path.name}',
                    ha='center', va='center', fontsize=11, color='gray',
                    transform=ax.transAxes)
            ax.set_title(title, fontsize=10, fontweight='bold', color='gray')
        ax.set_xticks([])
        ax.set_yticks([])
    
    save_path = REPO_DIR / "reports" / "barplot_aptr_density_left_newmask.png"
    fig.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.show()

# %% [markdown] QC Summary & Final Screenshot
# ---
# ## QC Summary
#
# | QC  | Check | Result |
# |-----|-------|--------|
# | 1c  | Brain mask covers entire brain (atlas-based), no holes | ☐ PASS / ☐ FAIL |
# | 2   | DWI-T2w registration — b0 brain edges match T2w in all 3 planes | ☐ PASS / ☐ FAIL |
# | 4   | FOD quality (directional glyphs in WM, small in ventricles) | ☐ PASS / ☐ FAIL |
# | 5   | ROIs in FOD space (aligned on FOD grid) | ☐ PASS / ☐ FAIL |
# | 6   | Segmentation (no overlaps, tiles seed ROI, SEG_MODE=9 for AP/APTR) | ☐ PASS / ☐ FAIL |
# | 7   | Whole-brain tractogram (fills brain volume) | ☐ PASS / ☐ FAIL |
# | 8   | Passthrough tracts (connect correct regions) | ☐ PASS / ☐ FAIL |
# | 9   | Cleaned tracts (tight bundles, outliers removed) | ☐ PASS / ☐ FAIL |
# | 10  | Tract metrics report generated (CSV + HTML) | ☐ PASS / ☐ FAIL |
# | 11  | AFD voxel map reasonable; ROI fixel metrics consistent | ☐ PASS / ☐ FAIL |
# | 11f | Fixel square ROI directions cluster correctly | ☐ PASS / ☐ FAIL |
#
# **Processed by:** _Your name_
# **Date:** _Date_
# **Notes:** _Any issues or observations_
