# MBRAIN Fixel-Based Analysis (FBA) - Task Tracker

## Date Started: 2026-02-08

---

## Current State Assessment

### What EXISTS:
| Item | Status | Details |
|------|--------|---------|
| Per-subject FODs (`wm_fod_norm.mif`) | DONE (12/12) | All subjects have individual FODs |
| Per-subject response functions | DONE (12/12) | `response_wm.txt`, `response_gm.txt`, `response_csf.txt` |
| Group-averaged response functions | DONE (12/12) | In `group_responses/`. Computed from ALL 12 subjects (verified from command_history) |
| Group-calibrated FODs (`wm_fod_norm_group.mif`) | **8/12** | AP/APTR done, **CTRL missing** |
| Per-subject preprocessed DWI | DONE (12/12) | AP/APTR: `dwi_flip_gradcheck_denoise_bias_anat.mif`, CTRL: both `dwi_preproc.mif` and `dwi_flip_gradcheck_denoise_bias_anat.mif` exist |
| Subject-space wholebrain tractograms | DONE (12/12) | `wholebrain/wholebrain_1000k.tck` + SIFT2 weights |
| Brain masks | DONE (12/12) | `mask.mif` for all subjects |
| Population FOD templates | NOT STARTED | No `templates/` directory exists |
| Subject-to-template warps | NOT STARTED | No warp files exist |
| Template fixel mask | NOT STARTED | |
| Fixel metrics (FD, FC, FDC) | NOT STARTED | |
| Template-space tractograms | NOT STARTED | |
| Statistical analysis | NOT STARTED | |

### CTRL _anat preprocessing bug (FIXED 2026-02-08):
- **Problem**: `dwi_flip_gradcheck_denoise_bias_anat.mif` had corrupted b-values for all 4 CTRL subjects
  - sub-129-921-R1F: 3 shells (114, 6652, 9330) instead of 2 shells (16, 1019)
  - sub-129-921-R1L1F: 4 shells (170, 7274, 10898, 11377)
  - sub-129-1014-L1M: 3 shells (93, 1049, 1477)
  - sub-138-1010-L1F: 1 shell (471)
- **Root cause**: The CTRL `dwi2anat_0GenericAffine.mat` (Jan 29) was NOT produced by 
  `antsRegistrationSyN.sh -d 3 -t r` (rigid). It contained a scaling/affine transform 
  (diagonal 0.34, 0.54, 0.67) instead of rigid (diagonal ~1.0). The `mrtransform -linear` 
  step then applied this scaling to gradient directions, corrupting b-values.
  Evidence: CTRL had no `Warped.nii.gz`/`InverseWarped.nii.gz` (standard `antsRegistrationSyN.sh` outputs).
  AP/APTR had correct near-identity rigid transforms.
- **Fix**: Re-ran `antsRegistrationSyN.sh -d 3 -t r` for all 4 CTRL subjects (steps 4j-4l from batch_process_flip1.py):
  1. ANTs rigid registration (b0_bias_brain.nii.gz → t2w_highres_flip.nii.gz)
  2. ConvertTransformFile + transformconvert (ITK → MRtrix format)
  3. mrtransform with corrected rigid transform
- **Result**: All 12 subjects now have consistent 2-shell structure: b≈16 (5 vols) + b≈1019 (40 vols)
- **BIDS source verification**: All 12 subjects have identical acquisition protocol in BIDS .bval files

### CTRL group FODs still needed:
- Group response averaging (`responsemean`) was already done correctly for all 12 subjects
- CTRL still need `wm_fod_norm_group.mif` generated with group response (Step 1b)

### Existing Scripts:
- `functions/compute_group_response.py` - Steps 1-3: collect responses, average, re-estimate FODs
  - Input DWI hardcoded to: `dwi_flip_gradcheck_denoise_bias_anat.mif`
  - CTRL subjects now have this file, so script should work with local paths
- `functions/analyze_roi_fixels.py` - ROI-based fixel extraction (INCOMPLETE: placeholder fixel-to-ROI mapping)
- `functions/fod_direction_analysis.py` - FOD directionality analysis within graft ROIs

---

## Analysis Pipeline Plan

### Why NOT a standard FBA pipeline (population template):
- **CTRL mice have a complete isocortex** (normal cortex)
- **AP/APTR mice are apallial** — absent/severely reduced dorsal telencephalon (cortex)
- These are fundamentally different brain morphologies
- MRtrix FBA docs (Step 10) warn: "avoid patients with excessive abnormalities compared to the rest of the population"
- A cross-group `population_template` would produce meaningless registrations in cortical regions
- **Per-group templates** (config1.yaml line 16) could be useful later, but NOT cross-group

### What IS shared across groups:
- **Common response function** (`responsemean` from all 12 subjects) — valid because it uses single-fiber WM voxels, which exist in all subjects regardless of cortical morphology
- **Group-calibrated FODs** (`wm_fod_norm_group.mif`) — computed with the shared response, making FD values comparable across subjects/groups

### Phase 1: Prerequisites ✅ COMPLETE
- [x] Step 1a: Group response functions computed (all 12 subjects, `responsemean`)
- [x] Step 1b: Group-calibrated FODs generated for all 12 subjects (`wm_fod_norm_group.mif`)

### Phase 2: ROI-Based Subject-Space Fixel Analysis ← CURRENT
- [x] Step 2a: Run `fod2fixel` on each subject's `wm_fod_norm_group.mif` in subject space
  - Output: fixel directories with FD (apparent fibre density) and directions per voxel
- [x] Step 2b: Extract fixels within shared ROIs (thalamus, striatum, paleocortex, cerebellum)
  - These ROIs exist in all 3 groups (from Allen atlas for CTRL, manual for AP/APTR)
  - Both already resampled to FOD space
- [x] Step 2c: Compute directional FD metrics per ROI (FD_total, FD_ap, FD_dv, FD_lm)
  - Use fixel directions to decompose FD into anterior-posterior, dorso-ventral, lateral-medial
- [x] Step 2d: Also extract fixels within seed regions (isocortex for CTRL, graft for AP/APTR)
  - These are group-specific ROIs (not directly comparable across groups)
  - Compare graft segments of AP vs APTR
- **OUTPUT**: `reports/fixel_all_rois.csv` — 120 rows (12 subjects × 5 ROIs × 2 hemispheres)
  - Columns: subject, group, roi, hemisphere, n_fixels, n_voxels, fd_total, fd_ap, fd_dv, fd_lm, fractions, fd_mean_per_fixel
  - Generated on 2026-02-25 with `python functions/analyze_roi_fixels.py --all-rois`

### Phase 3: Tract-Level Connectivity Metrics
- [ ] Step 3a: Use existing wholebrain tractograms + SIFT2 weights
  - Already available: `wholebrain/wholebrain_1000k.tck` + SIFT2 per subject
- [ ] Step 3b: Extract connectivity metrics per ROI pair (streamline count, mean AFD along tracts)
  - seed→thalamus, seed→striatum, seed→paleocortex, seed→cerebellum per hemisphere
- [ ] Step 3c: Also use segment-level tracts (already generated with tractography_segments.py)

### Phase 4: Statistical Comparison
- [ ] Step 4a: Aggregate fixel metrics into CSV (per subject, per ROI, per hemisphere)
- [ ] Step 4b: Group comparisons: CTRL vs AP, CTRL vs APTR, AP vs APTR
  - Use non-parametric tests (Mann-Whitney U, n=4 per group)
  - Kruskal-Wallis for 3-group omnibus test
- [ ] Step 4c: Generate connectivity report with plots and p-values

### Phase 5: Visualization
- [ ] Step 5a: Fixel overlay on T2w per representative subject per group
- [ ] Step 5b: Group bar plots of FD metrics
- [ ] Step 5c: Generate summary report

---

## Key Design Decisions
1. **NO cross-group population template** (morphologies too different: normal cortex vs apallial)
2. **Common response function** includes ALL 12 subjects (valid for WM)
3. **ROI-based subject-space analysis** instead of template-space FBA
4. **Shared ROIs** for comparison: thalamus, striatum, paleocortex, cerebellum (exist in all groups)
5. **Statistical contrasts**: CTRL vs AP, CTRL vs APTR, AP vs APTR (non-parametric, n=4/group)

---

## Subjects Reference

### CTRL (4 subjects - wild-type, normal cortex):
| Subject | Session |
|---------|---------|
| sub-129-921-R1F | ses-260115121231 |
| sub-129-921-R1L1F | ses-260116085433 |
| sub-129-1014-L1M | ses-260115113814 |
| sub-138-1010-L1F | ses-260115101441 |

### AP (4 subjects - apallial mutation, no transplant):
| Subject | Session |
|---------|---------|
| sub-596-C-OO-NmM | ses-240911103349 |
| sub-769-C-O5-R1M | ses-240911115733 |
| sub-769-C-OO-NmM | ses-240911111853 |
| sub-875-C-OO-NM | ses-240911100253 |

### APTR (4 subjects - apallial + cortical transplant):
| Subject | Session |
|---------|---------|
| sub-AP72-862-901-NmF | ses-251210093342 |
| sub-AP72-862-901-R1M | ses-251210100043 |
| sub-Apallial17NmF | ses-240816091633 |
| sub-Apallial20R1F | ses-240909094232 |

---

## Hemisphere Correction Log (2026-02-08)
- CTRL: Atlas ROIs had swapped L/R labels → FIXED (173 file pairs)
- AP/APTR: Had correct labels originally → confirmed no changes needed

## CTRL _anat Preprocessing Fix Log (2026-02-08)
- All 4 CTRL subjects had corrupted gradient tables in `dwi_flip_gradcheck_denoise_bias_anat.mif`
- Re-ran ANTs rigid registration (`antsRegistrationSyN.sh -d 3 -t r`) + transform conversion + mrtransform
- All 12 subjects now verified: 2 shells (b≈16, b≈1019), 45 volumes each
