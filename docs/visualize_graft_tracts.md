# Visualize Graft Tractography Results

Commands to visualize tractography from graft to target regions using mrview.
All commands use T2w highres as the background image.

---

## Setup

```bash
cd /Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc
```

---

## sub-596-C-OO-NmM

```bash
# All tracts together
mrview sub-596-C-OO-NmM/ses-240911103349/dwi/t2w_highres_flip.mif \
  -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Striatum_left_sift50.tck \
  -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Striatum_right_sift50.tck \
  -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Thalamus_left_sift50.tck \
  -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Thalamus_right_sift50.tck \
  -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Palecortex_left_sift50.tck \
  -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Palecortex_right_sift50.tck \
  -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_cerebellum_left_sift50.tck \
  -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_cerebellum_right_sift50.tck
```

```bash
# Striatum only
mrview sub-596-C-OO-NmM/ses-240911103349/dwi/t2w_highres_flip.mif \
  -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Striatum_left_sift50.tck \
  -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Striatum_right_sift50.tck
```

```bash
# Thalamus only
mrview sub-596-C-OO-NmM/ses-240911103349/dwi/t2w_highres_flip.mif \
  -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Thalamus_left_sift50.tck \
  -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Thalamus_right_sift50.tck
```

```bash
# Palecortex only
mrview sub-596-C-OO-NmM/ses-240911103349/dwi/t2w_highres_flip.mif \
  -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Palecortex_left_sift50.tck \
  -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Palecortex_right_sift50.tck
```

```bash
# Cerebellum only
mrview sub-596-C-OO-NmM/ses-240911103349/dwi/t2w_highres_flip.mif \
  -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_cerebellum_left_sift50.tck \
  -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_cerebellum_right_sift50.tck
```

---

## sub-Apallial17NmF

```bash
# All tracts together
mrview sub-Apallial17NmF/ses-240816091633/dwi/t2w_highres_flip.mif \
  -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Striatum_left_sift50.tck \
  -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Striatum_right_sift50.tck \
  -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Thalamus_left_sift50.tck \
  -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Thalamus_right_sift50.tck \
  -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Palecortex_left_sift50.tck \
  -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Palecortex_right_sift50.tck \
  -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_cerebellum_left_sift50.tck \
  -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_cerebellum_right_sift50.tck
```

```bash
# Striatum only
mrview sub-Apallial17NmF/ses-240816091633/dwi/t2w_highres_flip.mif \
  -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Striatum_left_sift50.tck \
  -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Striatum_right_sift50.tck
```

```bash
# Thalamus only
mrview sub-Apallial17NmF/ses-240816091633/dwi/t2w_highres_flip.mif \
  -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Thalamus_left_sift50.tck \
  -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Thalamus_right_sift50.tck
```

```bash
# Palecortex only
mrview sub-Apallial17NmF/ses-240816091633/dwi/t2w_highres_flip.mif \
  -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Palecortex_left_sift50.tck \
  -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Palecortex_right_sift50.tck
```

```bash
# Cerebellum only
mrview sub-Apallial17NmF/ses-240816091633/dwi/t2w_highres_flip.mif \
  -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_cerebellum_left_sift50.tck \
  -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_cerebellum_right_sift50.tck
```

---

## sub-769-C-OO-NmM

```bash
# All tracts together
mrview sub-769-C-OO-NmM/ses-240911111853/dwi/t2w_highres_flip.mif \
  -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Striatum_left_sift50.tck \
  -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Striatum_right_sift50.tck \
  -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Thalamus_left_sift50.tck \
  -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Thalamus_right_sift50.tck \
  -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Palecortex_left_sift50.tck \
  -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Palecortex_right_sift50.tck \
  -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_cerebellum_left_sift50.tck \
  -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_cerebellum_right_sift50.tck
```

```bash
# Striatum only
mrview sub-769-C-OO-NmM/ses-240911111853/dwi/t2w_highres_flip.mif \
  -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Striatum_left_sift50.tck \
  -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Striatum_right_sift50.tck
```

```bash
# Thalamus only
mrview sub-769-C-OO-NmM/ses-240911111853/dwi/t2w_highres_flip.mif \
  -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Thalamus_left_sift50.tck \
  -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Thalamus_right_sift50.tck
```

```bash
# Palecortex only
mrview sub-769-C-OO-NmM/ses-240911111853/dwi/t2w_highres_flip.mif \
  -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Palecortex_left_sift50.tck \
  -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Palecortex_right_sift50.tck
```

```bash
# Cerebellum only
mrview sub-769-C-OO-NmM/ses-240911111853/dwi/t2w_highres_flip.mif \
  -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_cerebellum_left_sift50.tck \
  -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_cerebellum_right_sift50.tck
```

---

## sub-AP72-862-901-R1M

```bash
# All tracts together
mrview sub-AP72-862-901-R1M/ses-251210100043/dwi/t2w_highres_flip.mif \
  -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Striatum_left_sift50.tck \
  -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Striatum_right_sift50.tck \
  -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Thalamus_left_sift50.tck \
  -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Thalamus_right_sift50.tck \
  -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Palecortex_left_sift50.tck \
  -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Palecortex_right_sift50.tck \
  -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_cerebellum_left_sift50.tck \
  -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_cerebellum_right_sift50.tck
```

```bash
# Striatum only
mrview sub-AP72-862-901-R1M/ses-251210100043/dwi/t2w_highres_flip.mif \
  -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Striatum_left_sift50.tck \
  -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Striatum_right_sift50.tck
```

```bash
# Thalamus only
mrview sub-AP72-862-901-R1M/ses-251210100043/dwi/t2w_highres_flip.mif \
  -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Thalamus_left_sift50.tck \
  -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Thalamus_right_sift50.tck
```

```bash
# Palecortex only
mrview sub-AP72-862-901-R1M/ses-251210100043/dwi/t2w_highres_flip.mif \
  -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Palecortex_left_sift50.tck \
  -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Palecortex_right_sift50.tck
```

```bash
# Cerebellum only
mrview sub-AP72-862-901-R1M/ses-251210100043/dwi/t2w_highres_flip.mif \
  -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_cerebellum_left_sift50.tck \
  -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_cerebellum_right_sift50.tck
```

---

## sub-AP72-862-901-NmF

```bash
# All tracts together
mrview sub-AP72-862-901-NmF/ses-251210093342/dwi/t2w_highres_flip.mif \
  -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Striatum_left_sift50.tck \
  -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Striatum_right_sift50.tck \
  -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Thalamus_left_sift50.tck \
  -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Thalamus_right_sift50.tck \
  -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Palecortex_left_sift50.tck \
  -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Palecortex_right_sift50.tck \
  -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_cerebellum_left_sift50.tck \
  -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_cerebellum_right_sift50.tck
```

```bash
# Striatum only
mrview sub-AP72-862-901-NmF/ses-251210093342/dwi/t2w_highres_flip.mif \
  -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Striatum_left_sift50.tck \
  -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Striatum_right_sift50.tck
```

```bash
# Thalamus only
mrview sub-AP72-862-901-NmF/ses-251210093342/dwi/t2w_highres_flip.mif \
  -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Thalamus_left_sift50.tck \
  -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Thalamus_right_sift50.tck
```

```bash
# Palecortex only
mrview sub-AP72-862-901-NmF/ses-251210093342/dwi/t2w_highres_flip.mif \
  -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Palecortex_left_sift50.tck \
  -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Palecortex_right_sift50.tck
```

```bash
# Cerebellum only
mrview sub-AP72-862-901-NmF/ses-251210093342/dwi/t2w_highres_flip.mif \
  -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_cerebellum_left_sift50.tck \
  -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_cerebellum_right_sift50.tck
```

---

## sub-Apallial20R1F

> **Note:** Cerebellum tracts had 0 tracks for this subject (small ROI size)

```bash
# All available tracts (no cerebellum)
mrview sub-Apallial20R1F/ses-240909094232/dwi/t2w_highres_flip.mif \
  -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Striatum_left_sift50.tck \
  -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Striatum_right_sift50.tck \
  -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Thalamus_left_sift50.tck \
  -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Thalamus_right_sift50.tck \
  -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Palecortex_left_sift50.tck \
  -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Palecortex_right_sift50.tck
```

```bash
# Striatum only
mrview sub-Apallial20R1F/ses-240909094232/dwi/t2w_highres_flip.mif \
  -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Striatum_left_sift50.tck \
  -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Striatum_right_sift50.tck
```

```bash
# Thalamus only
mrview sub-Apallial20R1F/ses-240909094232/dwi/t2w_highres_flip.mif \
  -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Thalamus_left_sift50.tck \
  -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Thalamus_right_sift50.tck
```

```bash
# Palecortex only
mrview sub-Apallial20R1F/ses-240909094232/dwi/t2w_highres_flip.mif \
  -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Palecortex_left_sift50.tck \
  -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Palecortex_right_sift50.tck
```

---

## sub-769-C-O5-R1M

```bash
# All tracts together
mrview sub-769-C-O5-R1M/ses-240911115733/dwi/t2w_highres_flip.mif \
  -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Striatum_left_sift50.tck \
  -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Striatum_right_sift50.tck \
  -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Thalamus_left_sift50.tck \
  -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Thalamus_right_sift50.tck \
  -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Palecortex_left_sift50.tck \
  -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Palecortex_right_sift50.tck \
  -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_cerebellum_left_sift50.tck \
  -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_cerebellum_right_sift50.tck
```

```bash
# Striatum only
mrview sub-769-C-O5-R1M/ses-240911115733/dwi/t2w_highres_flip.mif \
  -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Striatum_left_sift50.tck \
  -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Striatum_right_sift50.tck
```

```bash
# Thalamus only
mrview sub-769-C-O5-R1M/ses-240911115733/dwi/t2w_highres_flip.mif \
  -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Thalamus_left_sift50.tck \
  -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Thalamus_right_sift50.tck
```

```bash
# Palecortex only
mrview sub-769-C-O5-R1M/ses-240911115733/dwi/t2w_highres_flip.mif \
  -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Palecortex_left_sift50.tck \
  -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Palecortex_right_sift50.tck
```

```bash
# Cerebellum only
mrview sub-769-C-O5-R1M/ses-240911115733/dwi/t2w_highres_flip.mif \
  -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_cerebellum_left_sift50.tck \
  -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_cerebellum_right_sift50.tck
```

---

## Quick View All Subjects (All Tracts)

Run these one at a time to quickly check each subject:

```bash
cd /Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc

# 1. sub-596-C-OO-NmM
mrview sub-596-C-OO-NmM/ses-240911103349/dwi/t2w_highres_flip.mif -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Striatum_left_sift50.tck -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Striatum_right_sift50.tck -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Thalamus_left_sift50.tck -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Thalamus_right_sift50.tck -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Palecortex_left_sift50.tck -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Palecortex_right_sift50.tck -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_cerebellum_left_sift50.tck -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_cerebellum_right_sift50.tck

# 2. sub-Apallial17NmF
mrview sub-Apallial17NmF/ses-240816091633/dwi/t2w_highres_flip.mif -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Striatum_left_sift50.tck -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Striatum_right_sift50.tck -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Thalamus_left_sift50.tck -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Thalamus_right_sift50.tck -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Palecortex_left_sift50.tck -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_Palecortex_right_sift50.tck -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_cerebellum_left_sift50.tck -tractography.load sub-Apallial17NmF/ses-240816091633/dwi/tracks_Graft_cerebellum_right_sift50.tck

# 3. sub-769-C-OO-NmM
mrview sub-769-C-OO-NmM/ses-240911111853/dwi/t2w_highres_flip.mif -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Striatum_left_sift50.tck -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Striatum_right_sift50.tck -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Thalamus_left_sift50.tck -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Thalamus_right_sift50.tck -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Palecortex_left_sift50.tck -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_Palecortex_right_sift50.tck -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_cerebellum_left_sift50.tck -tractography.load sub-769-C-OO-NmM/ses-240911111853/dwi/tracks_Graft_cerebellum_right_sift50.tck

# 4. sub-AP72-862-901-R1M
mrview sub-AP72-862-901-R1M/ses-251210100043/dwi/t2w_highres_flip.mif -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Striatum_left_sift50.tck -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Striatum_right_sift50.tck -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Thalamus_left_sift50.tck -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Thalamus_right_sift50.tck -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Palecortex_left_sift50.tck -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_Palecortex_right_sift50.tck -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_cerebellum_left_sift50.tck -tractography.load sub-AP72-862-901-R1M/ses-251210100043/dwi/tracks_Graft_cerebellum_right_sift50.tck

# 5. sub-AP72-862-901-NmF
mrview sub-AP72-862-901-NmF/ses-251210093342/dwi/t2w_highres_flip.mif -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Striatum_left_sift50.tck -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Striatum_right_sift50.tck -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Thalamus_left_sift50.tck -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Thalamus_right_sift50.tck -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Palecortex_left_sift50.tck -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_Palecortex_right_sift50.tck -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_cerebellum_left_sift50.tck -tractography.load sub-AP72-862-901-NmF/ses-251210093342/dwi/tracks_Graft_cerebellum_right_sift50.tck

# 6. sub-Apallial20R1F (no cerebellum tracts)
mrview sub-Apallial20R1F/ses-240909094232/dwi/t2w_highres_flip.mif -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Striatum_left_sift50.tck -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Striatum_right_sift50.tck -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Thalamus_left_sift50.tck -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Thalamus_right_sift50.tck -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Palecortex_left_sift50.tck -tractography.load sub-Apallial20R1F/ses-240909094232/dwi/tracks_Graft_Palecortex_right_sift50.tck

# 7. sub-769-C-O5-R1M
mrview sub-769-C-O5-R1M/ses-240911115733/dwi/t2w_highres_flip.mif -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Striatum_left_sift50.tck -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Striatum_right_sift50.tck -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Thalamus_left_sift50.tck -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Thalamus_right_sift50.tck -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Palecortex_left_sift50.tck -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_Palecortex_right_sift50.tck -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_cerebellum_left_sift50.tck -tractography.load sub-769-C-O5-R1M/ses-240911115733/dwi/tracks_Graft_cerebellum_right_sift50.tck
```

---

## View with ROI Overlays

To also see the ROIs as overlays:

```bash
# Example for sub-596-C-OO-NmM with Graft and Striatum ROIs
mrview sub-596-C-OO-NmM/ses-240911103349/dwi/t2w_highres_flip.mif \
  -overlay.load sub-596-C-OO-NmM/ses-240911103349/manual_roi/dwi/roi_Graft_left_dwi.mif \
  -overlay.opacity 0.3 \
  -overlay.load sub-596-C-OO-NmM/ses-240911103349/manual_roi/dwi/roi_Striatum_left_dwi.mif \
  -overlay.opacity 0.3 \
  -tractography.load sub-596-C-OO-NmM/ses-240911103349/dwi/tracks_Graft_Striatum_left_sift50.tck
```
