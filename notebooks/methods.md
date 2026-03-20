# Methods

## Animals and experimental groups

Three groups of mice were used in this study (N = 4 per group, 12 total):
(1) **CTRL** — wild-type control mice with a fully developed, normal isocortex;
(2) **AP** — apallial mice carrying the *Emx1-Cre;Lhx2^{loxP/loxP}* mutation, which causes absence or severe reduction of the dorsal telencephalon, without any surgical intervention;
(3) **APTR** — mice carrying the same *Emx1-Cre;Lhx2^{loxP/loxP}* mutation that received an embryonic cortical tissue transplant. APTR animals were scanned approximately three months after transplantation.

## MRI acquisition

All imaging was performed ex vivo on a 7 T Bruker scanner (Stanford 70/40, Stanford University) running ParaVision 360.3.5 software and equipped with an RF CP F2.5 300 1 coil. Medical air was used as the breathing gas to minimize the effects of high O₂ concentration on relaxation.

**Diffusion-weighted imaging (DWI).** Data were acquired with a 4-segment diffusion-weighted EPI sequence (Bruker:DtiEpi) using the following parameters: b = 1000 s/mm², 40 non-collinear diffusion-encoding directions plus 5 b = 0 images (45 volumes total), diffusion timing δ/Δ = 4.5/10.3 ms, TR/TE = 3500/21 ms, in-plane resolution 150 × 150 µm², slice thickness 500 µm, FOV = 18 × 15 mm², matrix = 120 × 100, bandwidth = 250 kHz, second-order B₀ map shimming, three FOV saturation bands.

**T2-weighted anatomical imaging.** Two RARE (Rapid Acquisition with Relaxation Enhancement) sequences were acquired per animal: (i) a low-resolution scan (200 × 200 × 500 µm³) for co-registration reference, and (ii) a high-resolution scan (100 × 100 × 300 µm³; RARE factor = 4, TR/TE = 4000/30 ms, 2 averages, flip angle = 180°) used as the primary anatomical reference for visualization and registration.

## Image preprocessing

All diffusion and anatomical data were preprocessed using MRtrix3 (v3.0.3; Tournier et al., 2019) and ANTs (Avants et al., 2009). The processing pipeline was applied identically to all 12 subjects.

### Orientation standardization

Raw images were reoriented to match the Allen Common Coordinate Framework (Allen CCF) via a two-step transformation ("flip1"): (i) axis permutation (axes 2, 1, 0) to map scanner axes to anterior–posterior (AP), dorsal–ventral (DV), and medial–lateral (ML) axes; and (ii) a flip along axis 1 (Y) to align the dorsal–ventral direction with the atlas convention. The same transformation was applied to the 4-D DWI data (axes 2, 1, 0, 3) and the high-resolution T2w image.

### DWI preprocessing

DWI preprocessing proceeded through the following steps:
1. **Gradient verification.** The diffusion gradient table was checked and corrected for orientation consistency using `dwigradcheck`.
2. **Denoising.** Thermal noise was removed using MP-PCA denoising (Veraart et al., 2016) with a 5 × 5 × 5 voxel kernel and the Exp2 estimator.
3. **Brain mask estimation.** An initial brain mask was created from the denoised DWI data.
4. **Bias field correction.** B₁ field inhomogeneity was corrected using the N4 algorithm (Tustison et al., 2010) from ANTs (B-spline fitting distance = 200, convergence [10, 1e-6], sigma = 4).
5. **DWI-to-T2w registration.** The bias-corrected mean b = 0 image was rigidly registered to the high-resolution T2w anatomical using ANTs (`antsRegistrationSyN.sh -t r`). The resulting affine transform was converted to MRtrix format and applied to the full 4-D DWI volume, bringing all diffusion data into T2w anatomical space.
6. **Final brain mask.** A refined brain mask was computed from the registered DWI data for use in all downstream analyses.

### Atlas-based brain mask generation

To ensure consistent brain masking across all subjects, a group-specific atlas-based masking approach was adopted. One representative subject per group was selected (CTRL: sub-129-921-R1F; AP: sub-596-C-OO-NmM; APTR: sub-AP72-862-901-NmF) and its brain mask was manually corrected using the mrview ROI editor. This edited mask was then propagated to the remaining subjects within the same group via ANTs SyN nonlinear registration: the atlas T2w image was registered to each target T2w image, and the atlas mask was warped using nearest-neighbor interpolation. The largest connected component was retained after morphological cleaning. Masks were resampled from T2w resolution to DWI resolution using nearest-neighbor interpolation prior to all downstream processing.

## Fiber orientation distribution estimation

Multi-tissue response functions were estimated per subject using the dhollander algorithm (Dhollander et al., 2016) for white matter (WM), grey matter (GM), and cerebrospinal fluid (CSF) compartments. For cross-subject comparability, response functions were averaged across all 12 subjects to create group-level response functions. FODs were then re-estimated per subject using multi-shell multi-tissue constrained spherical deconvolution (MSMT-CSD; Jeurissen et al., 2014) with the group-averaged WM and CSF response functions. Intensity normalization was performed using `mtnormalise` (Raffelt et al., 2017) with WM and CSF compartments only (three-tissue normalization including GM was unstable for mouse brain data). When normalization failed, unnormalized FODs were used as a fallback. The resulting group-normalized WM FOD images (`wm_fod_norm_group.mif`) were used for all subsequent analyses.

## Region of interest definition

### CTRL subjects — Allen CCF atlas ROIs

For CTRL subjects, regions of interest were derived from the Allen CCF atlas (Allen Institute for Brain Science). The Allen average template was nonlinearly registered to each subject's T2w image using ANTs SyN, and the warp was applied to the Allen annotation volume using nearest-neighbor interpolation. ROIs were extracted by collecting all descendant structure IDs under the following parent structures: isocortex (ID 315), piriform complex (PIR 961 + ENT 909 + PAA 788 + TR 566), striatum (STRd 485 + STRv 493), thalamus (ID 549), and cerebellum (ID 512). Each ROI was split into left and right hemispheres along the midline of axis 2.

### AP/APTR subjects — manual ROIs

For AP and APTR subjects, ROIs were manually segmented in 3D Slicer on the high-resolution T2w images. Five bilateral structures were delineated: Graft (the abnormal cortical tissue in AP or the transplanted tissue in APTR), Striatum, Thalamus, Palecortex, and cerebellum. Manual ROIs were transformed to flip1 orientation and resampled to FOD space using nearest-neighbor interpolation.

### Seed ROI segmentation

For cross-group comparison, the seed ROI (isocortex for CTRL, Graft for AP/APTR) was divided along the anterior–posterior axis into three equal segments (anterior, central, posterior), yielding bilateral dorsal sub-regions. For AP/APTR subjects, a 3 × 3 segmentation scheme (anterior/central/posterior × dorsal/medial/ventral = 9 segments) was also available; the dorsal segments were used for the primary analysis. For AP/APTR, eroded ("shrunk") versions of the seed segments were generated by applying binary erosion (1 voxel, 6-connectivity; SciPy) to reduce partial-volume effects at ROI boundaries.

## Whole-brain tractography

Whole-brain tractograms were generated using the iFOD2 probabilistic algorithm (Tournier et al., 2010) seeded from the brain mask. Tracking parameters were: step size = 0.02 mm, maximum angle between successive steps = 60°, FOD amplitude cutoff = 0.1, minimum streamline length = 2 mm, maximum streamline length = 15 mm. One million streamlines were generated per subject.

## SIFT2 and passthrough tract extraction

Streamline weights were estimated using SIFT2 (Smith et al., 2015) to achieve a biologically meaningful estimation of fiber connectivity strength. Passthrough tracts were extracted from the whole-brain tractograms using `tckedit` with `-include` flags for ipsilateral seed-segment-to-target ROI pairs. Targets included thalamus, striatum, piriform/olfactory cortex, and cerebellum (the latter requiring an additional thalamus waypoint to constrain cortico-thalamo-cerebellar pathways).

## Tract cleaning

Outlier streamlines were removed using the AFQ_removeFiberOutliers function (Yeatman et al., 2012) from the Automated Fiber Quantification (AFQ) toolbox, called via MATLAB (R2024a, MathWorks). Streamlines deviating more than 3 standard deviations from the tract core in either Gaussian distance or fiber length were iteratively removed (up to 5 iterations, 100 resampling nodes per fiber, mean core method). Tracts with fewer than 20 streamlines were excluded from cleaning.

## Fixel-based analysis

FODs were decomposed into discrete fixels using `fod2fixel`, yielding per-fixel apparent fiber density (AFD) and peak amplitude maps. The voxel-wise total AFD was computed by summing all fixel contributions via `fixel2voxel`. Directional fiber density within each ROI was decomposed into AP, DV, and ML components for ROI-level microstructural characterization.

### Slicewise fixel direction analysis

To characterize the spatial organization of fiber orientations within the seed ROI, a slicewise fixel analysis was performed. For each coronal slice through the seed ROI (isocortex for CTRL, Graft for AP/APTR), fixel peak direction vectors were extracted from the group-normalized FOD using `fod2fixel` followed by `fixel2peaks`. Fixel directions were clustered into two populations using K-means clustering (k = 2) applied to the absolute values of the direction vectors (to handle antipodal symmetry). For each cluster, the mean direction was computed using the scatter-matrix eigenvector method, which correctly averages axial data by finding the principal eigenvector of $\mathbf{S} = \sum_i \mathbf{v}_i \mathbf{v}_i^\top$. Cluster directions were scaled by the mean FOD peak amplitude to produce amplitude-weighted direction vectors.

### Fixel dominance analysis

The dominant fiber direction in each voxel was determined by the axis with the largest absolute component of the fixel direction vector (AP, DV, or ML). The percentage of fixels dominated by each direction was calculated per subject and averaged per group.

### Direction vector RMSD analysis

Within-group consistency of fiber orientations was quantified by the root-mean-square deviation (RMSD) of direction vectors across subjects, computed per slice and overall:

$$\text{RMSD} = \sqrt{\frac{1}{N}\sum_{i=1}^{N} \|\mathbf{v}_i - \bar{\mathbf{v}}\|^2}$$

where $\mathbf{v}_i$ is the K-means cluster 1 direction vector (scaled by FOD peak amplitude) for subject $i$, and $\bar{\mathbf{v}}$ is the group mean vector. The mean vector norm $\|\bar{\mathbf{v}}\|$ and the ratio RMSD/$\|\bar{\mathbf{v}}\|$ were also computed to contextualize variability relative to signal magnitude.

## Visualization

### Radar plots

Fixel direction vectors were visualized as radar (polar) plots with three axes corresponding to the AP, DV, and ML anatomical directions. Each radar plot displays one coronal slice through the seed ROI, with individual mice shown as separate traces using group-specific color schemes (grey shades for CTRL, magenta shades for AP, blue shades for APTR) and per-mouse marker symbols. A thin radial grid with three concentric rings provides the amplitude scale. The amplitude axis was unified across all groups and slices for direct comparison.

### Composite figures

Three-group comparison figures were generated showing fixel direction radar plots arranged in a grid with groups as rows and coronal slices as columns, enabling visual comparison of directional fiber organization across experimental conditions and along the anterior–posterior axis.

## Software

All analyses were implemented in Python (≥3.10) using NumPy, SciPy, scikit-learn, nibabel, and Matplotlib. DWI processing was performed with MRtrix3 v3.0.3 (Tournier et al., 2019). Image registration used ANTs (Avants et al., 2009). Atlas ROI extraction used the Allen CCF (Allen Institute for Brain Science). Manual ROI segmentation was performed in 3D Slicer. Tract cleaning used MATLAB R2024a with the AFQ toolbox (Yeatman et al., 2012). BIDS-formatted data were converted from Bruker ParaVision format using dcm2niix (Li et al., 2016).

## References

- Avants, B. B., Tustison, N. J., Song, G., Cook, P. A., Klein, A., & Gee, J. C. (2009). A reproducible evaluation of ANTs similarity metric performance in brain image registration. *NeuroImage*, 54(3), 2033–2044.
- Dhollander, T., Raffelt, D., & Connelly, A. (2016). Unsupervised 3-tissue response function estimation from single-shell or multi-shell diffusion MR data without a co-registered T1 image. *ISMRM Workshop on Breaking the Barriers of Diffusion MRI*, 5.
- Jeurissen, B., Tournier, J.-D., Dhollander, T., Connelly, A., & Sijbers, J. (2014). Multi-tissue constrained spherical deconvolution for improved analysis of multi-shell diffusion MRI data. *NeuroImage*, 103, 411–426.
- Li, X., Morgan, P. S., Ashburner, J., Smith, J., & Rorden, C. (2016). The first step for neuroimaging data analysis: DICOM to NIfTI conversion. *Journal of Neuroscience Methods*, 264, 47–56.
- Raffelt, D. A., Tournier, J.-D., Smith, R. E., Vaughan, D. N., Jackson, G., Ridgway, G. R., & Connelly, A. (2017). Investigating white matter fibre density and morphology using fixel-based analysis. *NeuroImage*, 144, 58–73.
- Smith, R. E., Tournier, J.-D., Calamante, F., & Connelly, A. (2015). SIFT2: Enabling dense quantitative assessment of brain white matter connectivity using streamlines tractography. *NeuroImage*, 119, 338–351.
- Tournier, J.-D., Calamante, F., & Connelly, A. (2010). Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. *Proceedings of the International Society for Magnetic Resonance in Medicine*, 18, 1670.
- Tournier, J.-D., Smith, R. E., Raffelt, D. A., Tabbara, R., Dhollander, T., Pietsch, M., ... & Connelly, A. (2019). MRtrix3: A fast, flexible and open software framework for medical image processing and visualisation. *NeuroImage*, 202, 116137.
- Tustison, N. J., Avants, B. B., Cook, P. A., Zheng, Y., Egan, A., Yushkevich, P. A., & Gee, J. C. (2010). N4ITK: improved N3 bias correction. *IEEE Transactions on Medical Imaging*, 29(6), 1310–1320.
- Veraart, J., Novikov, D. S., Christiaens, D., Ades-Aron, B., Sijbers, J., & Fieremans, E. (2016). Denoising of diffusion MRI using random matrix theory. *NeuroImage*, 142, 394–406.
- Yeatman, J. D., Dougherty, R. F., Myall, N. J., Wandell, B. A., & Feldman, H. M. (2012). Tract profiles of white matter properties: automating fiber-tract quantification. *PLoS ONE*, 7(11), e49790.
