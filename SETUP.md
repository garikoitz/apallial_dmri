# Quick Setup Guide for the apallial_dmri Pipeline

## Quick Start (New Computer)

### Step 1: Clone and create environment

```bash
git clone https://github.com/garikoitz/apallial_dmri.git
cd apallial_dmri
micromamba env create -f mbrain_env.yml
micromamba activate mbrain
```

### Step 2: Download the data

```bash
python utils/setup_data.py
```

This downloads and extracts all 4 OSF components to `local/` and `data/`.
**Estimated time:** 20-40 minutes (depends on internet speed).

To download only specific components:

```bash
python utils/setup_data.py --components local          # atlas + ROIs only
python utils/setup_data.py --components pre_bids       # BIDS dataset only
python utils/setup_data.py --components pre_preproc    # pre-mask derivatives only
python utils/setup_data.py --components post           # post-mask derivatives only
python utils/setup_data.py --check-only                # check what's missing
```

### Step 3: Run the notebook

The pipeline notebook is a jupytext percent-format `.py` file
(`notebooks/pipeline_groups_new_mask.py`). The recommended way to run it is
**inside VS Code as an interactive notebook** (no `jupyter notebook` server
needed):

1. Install the **Jupyter** and **Python** VS Code extensions.
2. Select the `mbrain` interpreter: `Cmd+Shift+P` →
   `Python: Select Interpreter` → choose `mbrain`.
3. Open `notebooks/pipeline_groups_new_mask.py` in VS Code.
4. Click **"Open in Notebook"** at the top-right of the editor (or run the
   command `Jupyter: Open as Notebook`). Each `# %%` block becomes a runnable
   cell; run with `Shift+Enter`. Edits save back to the same `.py` file.

> `jupytext` is already in `mbrain_env.yml`, so the env created in Step 1 has
> everything needed.

Alternative (browser-based):

```bash
micromamba run -n mbrain jupyter lab notebooks/pipeline_groups_new_mask.py
```

## Basic Replication Workflow (Recommended)

### 1. Edit config1.yaml first

Update these fields for your machine and run scope:

- paths.repo_dir
- paths.data_base
- paths.matlab_path
- groups (subject/session roster)

### 2. Install/check data with setup_data.py

```bash
python utils/setup_data.py --check-only
python utils/setup_data.py --components local pre_bids pre_preproc post
```

### 3. Open the pipeline notebook in VS Code

Open notebooks/pipeline_groups_new_mask.py and use Open in Notebook with the mbrain kernel.

### 4. Run notebook steps by starting point

- If data/derivatives/mrtrix_newmask is already populated: skip Steps 1-4, start at Step 0.
- If starting from BIDS + local only: run Steps 1, 2, 3, 4, then run Step 0 and continue.

## What Happens Automatically

When you run `python utils/setup_data.py`:

1. Downloads the 4 OSF component ZIPs to `~/.mbrain_temp/`
2. Extracts them to the repository root (`local/` and `data/`)
3. Creates marker files (`.mbrain_component_*`) to track what's installed
4. Validates critical files are present
5. You're ready to run the notebook

**Note:** `local/` and `data/` are gitignored — they won't be committed.

## Three Paths Forward

### Path A: Fast-Track (Recommended for Testing)
**If you want to see results quickly — install all 4 components:**

```bash
python utils/setup_data.py --components local pre_bids pre_preproc post
```

Then in the notebook:
```
Step 1 (verify data)
    ↓
Skip Steps 2, 3, 4
    ↓
Step 0 onwards (new-mask pipeline)
    ↓
Run in ~1 hour
```

The `data/derivatives/mrtrix_newmask/` directory is already pre-computed in the `pre_preproc` + `post` components.

### Path B: Full Reproduction (Recommended for Methods)
**If you want to reprocess from BIDS — install local + pre_bids + pre_preproc:**

```bash
python utils/setup_data.py --components local pre_bids pre_preproc
```

Then in the notebook:
```
Step 1 (verify data)
    ↓
Step 2 (T2w flip + DWI preprocessing) — 30-60 min per subject
    ↓
Step 3 (FOD estimation + ROI resampling) — 30-45 min per subject
    ↓
Step 4 (Segmentation + tractography) — 2-4 hours per subject
    ↓
Step 0 onwards (new-mask pipeline)
    ↓
Run in 12-36 hours total (depending on machine)
```

Uses `data/bids/` to process everything fresh, then syncs outputs into `data/derivatives/mrtrix_newmask/` for Step 0+.

### Path C: Bring Your Own Data
**If you have your own subjects:**

1. Follow Path B to understand the preprocessing pipeline
2. Adapt Steps 1-4 for your subjects in data/bids/
3. Run Steps 1-4 to generate flip1 derivatives
4. Then run Steps 0+ for mask analysis

## Troubleshooting

### Q: Auto-download keeps failing
**A:** Try manual download:
```bash
python utils/setup_data.py --check-only    # See what's missing
python utils/setup_data.py --force         # Force re-download
```

### Q: Takes too long to download / Interrupted
**A:** Download is resumable. Just run again:
```bash
python utils/setup_data.py
# If it times out, re-run and it will pick up where it left off
```

### Q: Want to skip download and use pre-existing data
**A:** If you already have data elsewhere, create symlinks:
```bash
# If data is on external drive:
ln -s /path/to/external/mbrain_data/local ./local
ln -s /path/to/external/mbrain_data/data ./data
```

### Q: Not enough disk space
**A:** Check available space:
```bash
df -h
# You need ~15 GB (6 GB ZIP + 9 GB unzipped)
```

If low on space, extract ZIP to external drive and symlink from there.

### Q: Which Python/Conda environment?
**A:** Use the environment specified in `mbrain_env.yml`:
```bash
micromamba env create -f mbrain_env.yml
micromamba activate mbrain
micromamba run -n mbrain jupyter lab notebooks/pipeline_groups_new_mask.py
```

## Data Contents

The data is split into four ZIP components:

| Component | Installs to | Contents | Purpose |
|-----------|-------------|----------|---------|
| `local` | `local/` | Allen CCF atlas + manual ROIs | Required reference inputs |
| `pre_bids` | `data/bids/` | BIDS dataset (dataset_description.json, participants.tsv, etc.) | Raw data for full reproduction |
| `pre_preproc` | `data/derivatives/mrtrix_newmask/` | Pre-mask derivatives (DWI, masks, FODs, ROIs) | Start from mask creation |
| `post` | `data/derivatives/mrtrix_newmask/` | Post-mask derivatives (tractograms, SIFT2, cleaned tracts) | Skip downstream recomputation |

**Total size:** ~6 GB compressed, ~9 GB unzipped (varies by component selection)

## For Paper Reviewers / Reproducibility

To verify results, you have two options:

**Quick verification (1-2 hours):**
- Install all 4 components: `python utils/setup_data.py`
- Skip Steps 1-4 and use pre-computed outputs
- Run Step 0 onwards to generate new masks and verify outputs

**Full reproducibility (12-36 hours):**
- Install `local` + `pre_bids` + `pre_preproc`: `python utils/setup_data.py --components local pre_bids pre_preproc`
- Run Steps 1 through 4 to reproduce all preprocessing from BIDS
- Run Step 0 onwards to generate new masks
- Verify every intermediate output against published derivatives

## Reference

- **Main notebook:** `notebooks/pipeline_groups_new_mask.py`
- **Processing guide:** `README.md` (detailed methods section)
- **Data packaging:** `docs/OSF_DATA_PACKAGING.md`
- **Configuration:** `config1.yaml` (all pipeline parameters)

## Citation

If you use this pipeline/data, cite:

```bibtex
@article{glerma_2026_mbrain,
  title={MBRAIN: Multi-subject Mouse Brain Rodent Analysis In Neuroimaging},
  author={Glerma, Garikoitz and ...},
  journal={...},
  year={2026}
}

@misc{glerma_2026_testdata,
  title={MBRAIN Test Data (12 subjects, BIDS + derivatives)},
  author={Glerma, Garikoitz},
  year={2026},
  url={https://osf.io/ABC123/},
  publisher={Open Science Framework}
}
```

---

**Need help?** See the main [README.md](README.md) or open an issue on GitHub.
