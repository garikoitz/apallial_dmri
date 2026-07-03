# Quick Start: Create & Distribute Test Data

This guide walks through the next steps to create the data package and upload to OSF.io.

## Prerequisites

- Access to the full data directory (`/export/home/.../mbrain/data/bids/` + `/export/home/.../mbrain/data/derivatives/mrtrix_newmask/`)
- Linux machine with ~30 GB free disk space (for staging + ZIP)
- ~30-60 minutes (20 min package creation + 30 min OSF upload)


```bash
bash utils/validate_setup.sh
```

Should see: ✅ All checks passed!

## Step 1: Create Data Package ZIP

**On the machine with data access (e.g., lmc02 with `/export/home/` mount):**

```bash
 This guide walks through the next steps to create the three data packages and upload them to OSF.io.

# Edit the script to point to your data directory (if not /export/home/...)
nano utils/create_data_package.sh
# Update line 8: DATA_BASE="/export/home/.../mbrain"

# Create the ZIP (takes ~20 minutes)
bash utils/create_data_package.sh

# Verify output
ls -lh /tmp/mbrain_test_data.zip
cat /tmp/mbrain_test_data.zip.sha256
```

Note: the packaging script now packages `data/derivatives/mrtrix_newmask/`
directly, so manually edited masks already present there
(`mask_new_edited.mif`, and `.nii.gz` when present) are included as-is.

**Expected output:**
```
✓ Verifying required files...
✓ Creating ZIP (phase 3/7)...
✓ Validating ZIP integrity (phase 4/7)...
✓ Computing SHA256 (phase 5/7)...
-----
5.8 GB ZIP file created
SHA256: a1b2c3d4e5f6...
```

## Step 2: Upload to OSF.io

**On any machine with internet + web browser:**

Follow the exact checklist in:
```
docs/OSF_WORKFLOW_CHECKLIST.md
```
 3. Upload the three ZIPs to the project
Key steps:
 4. Get the three download URLs and store them as component-specific links
1. Create free OSF.io account
2. Create new project (name: `mbrain-test-data`)
3. Upload `/tmp/mbrain_test_data.zip` to project
4. Get download URL: `https://osf.io/ABC123/download`

**Time:** ~30 min (depending on upload speed)

## Step 3: Configure Repository

Update the OSF URL in the setup script:

```bash
# Edit setup_data.py
 # Set component URLs:
 # MBRAIN_URL_LOCAL=...
 # MBRAIN_URL_PRE=...
 # MBRAIN_URL_POST=...
nano utils/setup_data.py

# Change line 12:
# DEFAULT_OSF_URL = "https://osf.io/ABC123/download"
#                                     ^^^^^^
#                                     Use your project ID from Step 2
```

## Step 4: Test Configuration

Verify auto-download works:

```bash
# Dry run: check structure without downloading
python utils/setup_data.py --check-only

# You should see:
# Missing: local/atlas/
# Missing: local/manual_roi_...
# 
# (This means download is ready to trigger)

# Actual download (on machine WITHOUT local/ folder):
# python utils/setup_data.py
# 
# Watch for:
# ✓ Downloading from OSF...
# ✓ Extracting ZIP...
# ✓ Validating files...
```

## Step 5: Commit & Push to GitHub

- Add utils/create_data_package.sh for packaging data
- Add utils/validate_setup.sh for setup verification
- Add Pre-Steps (-2, -1a, -1b, -1c) to main analysis notebook
- Update documentation (SETUP.md, osf_*.md)
- Data: 5.8 GB archive, auto-downloads on first run
- Supports both fast-track and full reproducibility
"

git push origin main
```

## Step 6: Paper Integration

Add to your Methods section:

```
### Data Availability

Test data (N=12 subjects) are available on OSF.io [1]:
https://osf.io/ABC123/

Raw DWI is publicly available via Flywheel (see Supplementary Methods).
Preprocessed diffusion derivatives (flip1 processing) are ~5.8 GB and 
auto-download on first notebook run. Users may skip preprocessing by using 
the fast-track mode (1 hour) or re-run from BIDS data (12-36 hours).

All analysis code is available at: https://github.com/YOUR_REPO

[1] OSF DOI: https://doi.org/10.17605/OSF.IO/ABC123 (if you register the 
    project; see OSF Dashboard → Project Settings → Registration)
```

Then cite in references:
```
[1] Glerma M, et al. (2024). Mbrain analysis code and test data.

## Troubleshooting
- Edit line 8 in `utils/create_data_package.sh`
- Ensure `DATA_BASE` points to a directory that contains `data/` subdirectory
- Consider: `ssh lmc02` → run script → `scp` back to local

- ZIP can be split: `split -b 1G mbrain_test_data.zip` (upload parts separately)

- Or manually run: `python utils/setup_data.py`

### Auto-download URL not working
- Verify you updated `DEFAULT_OSF_URL` in line 12 of `setup_data.py`
- Test: `python utils/setup_data.py --check-only`
- If still broken, check URL format: `https://osf.io/ABC123/download`

## Reference Files

- **Quick reference:** This file
- **User setup:** `SETUP.md`
- **OSF workflow:** `docs/OSF_WORKFLOW_CHECKLIST.md`
- **Technical details:** `docs/OSF_DATA_PACKAGING.md`
- **Implementation summary:** `IMPLEMENTATION_SUMMARY.md`

## Questions?

Refer to the appropriate documentation:
- **"How do I use this?"** → `SETUP.md`
- **"What's the technical approach?"** → `docs/OSF_DATA_PACKAGING.md`
- **"What's the step-by-step OSF process?"** → `docs/OSF_WORKFLOW_CHECKLIST.md`
- **"What code changed?"** → `IMPLEMENTATION_SUMMARY.md`
