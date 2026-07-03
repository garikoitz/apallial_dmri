# Quick Setup Guide for MBRAIN Pipeline

## One-Line Start

```bash
git clone https://github.com/garikoitz/apallial_dmri.git
cd apallial_dmri
jupyter notebook notebooks/pipeline_groups_new_mask.py
# Pre-Step -2 will auto-download data on first run
```

## What Happens Automatically

When you run the notebook for the first time:

1. **Pre-Step -2** checks if `local/` directory exists
2. If missing, it automatically downloads the selected OSF component ZIP(s)
3. Extracts the ZIPs to your repository
4. Validates all critical files
5. You're ready to run Pre-Steps or skip to Step 0

**Estimated time:** 20-40 minutes (depends on internet speed)

## Three Paths Forward

### Path A: Fast-Track (Recommended for Testing)
**If you want to see results quickly:**

```
Pre-Step -2 (verify data)
    ↓
Skip Pre-Steps -1a, -1b, -1c
    ↓
Step 0 onwards (new-mask pipeline)
    ↓
Run in ~1 hour
```

The `data/derivatives/mrtrix_newmask/` directory is already pre-computed.

### Path B: Full Reproduction (Recommended for Methods)
**If you want to reprocess from BIDS:**

```
Pre-Step -2 (verify data)
    ↓
Pre-Step -1a (T2w flip + DWI preprocessing) — 30-60 min per subject
    ↓
Pre-Step -1b (FOD estimation + ROI resampling) — 30-45 min per subject
    ↓
Pre-Step -1c (Segmentation + tractography) — 2-4 hours per subject
    ↓
Step 0 onwards (new-mask pipeline)
    ↓
Run in 12-36 hours total (depending on machine)
```

Uses `data/bids/` to process everything fresh, then syncs outputs into `data/derivatives/mrtrix_newmask/` for Step 0+.

### Path C: Bring Your Own Data
**If you have your own subjects:**

1. Follow Path B to understand the preprocessing pipeline
2. Adapt Pre-Steps for your subjects in `data/bids/`
3. Run Pre-Steps to generate flip1 derivatives
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
micromamba create -f mbrain_env.yml
micromamba activate mbrain
jupyter notebook notebooks/pipeline_groups_new_mask.py
```

## Data Contents

The data is split into three ZIPs:

| Path | Contents | Purpose |
|------|----------|---------|
| `local` | Allen CCF atlas + manual ROIs | Required reference inputs |
| `pre` | BIDS + pre-mask mrtrix_newmask files | Start from BIDS or skip to mask creation |
| `post` | Remaining mrtrix_newmask outputs | Skip downstream recomputation |

**Total size:** split across three files; the exact size depends on the selected component(s)

## For Paper Reviewers / Reproducibility

To verify results, you have two options:

**Quick verification (1-2 hours):**
- Install `local` + `pre` + `post`, then skip Pre-Steps and use pre-computed outputs
- Run Step 0 onwards to generate new masks and verify outputs

**Full reproducibility (12-36 hours):**
- Install `local` + `pre`, run Pre-Steps -2 through -1c to reproduce all preprocessing from BIDS
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
