# RamanPCACluster

Unsupervised PCA + k-means clustering of Raman spectra (`.dpt` format),
as a MATLAB script and an interactive App.

## Requirements

- MATLAB with **Statistics and Machine Learning Toolbox** (`pca`, `kmeans`,
  `silhouette`)
- **Signal Processing Toolbox** (`sgolayfilt`)
- No external dependencies beyond MATLAB itself: `backcor.m`, `airPLS.m`,
  `snip.m`, `apls.m`, and `readdpt.m` are included in this folder.

## Files

- `PCAClusterApp.m` — interactive app: pick a folder of `.dpt` spectra,
  set the spectral range, preprocessing, and number of clusters, run the
  analysis, inspect the results across tabs (including a hierarchical
  clustering dendrogram, and optional 80/85/90% confidence ellipses on
  the PCA scatter plots), and export them.
- `PCA_kmeans_analysis.m` — the same pipeline as a plain script (edit the
  top of the file to point at a different `Spectra/` folder), useful for
  batch/reproducible runs outside the GUI.
- `loadRamanSpectra.m` — loads every `.dpt` file in a folder onto a common
  wavenumber grid; also derives a group label per spectrum from its
  filename (e.g. `CT_P_11b.0.dpt` → `CT_P`), used only for
  validating/plotting the clustering, never for the clustering itself.
- `preprocessSpectra.m` — baseline removal (choice of backcor / airPLS /
  SNIP / APLS, optional), Savitzky-Golay smoothing (optional), and
  normalization (area / max / SNV / none) -- the same four baseline
  methods RamanFitApp offers.
- `adjustedRandIndex.m` — Adjusted Rand Index between two partitions
  (e.g. k-means clusters vs. filename-derived groups), corrected for
  chance agreement.
- `backcor.m`, `airPLS.m`, `snip.m`, `apls.m`, `readdpt.m` —
  third-party/personal helper functions, copied in for self-containment
  (see below).

## Pipeline

1. **Load**: every `.dpt` spectrum in a folder, stacked into one matrix
   (all files must share the same wavenumber grid).
2. **Restrict** to a chosen spectral range (optional).
3. **Preprocess**: baseline removal (airPLS) + light smoothing +
   normalization, so PCA isn't dominated by fluorescence background or
   overall signal-intensity differences between measurements.
4. **PCA**: scree plot; enough PCs to explain ≥95% of variance (capped at
   10) are kept for clustering.
5. **Choose k**: elbow (within-cluster sum of squares) and silhouette
   analysis over k = 2..8.
6. **k-means** on the retained PCA scores, with a fixed or
   silhouette-suggested k.
7. **Compare** the resulting clusters against each spectrum's
   filename-derived group (contingency table + Adjusted Rand Index) --
   informational only, not used to guide the clustering.
8. **Dendrogram** (App only): Ward-linkage hierarchical clustering on the
   same retained PCA scores, as an alternative view of cluster structure
   that doesn't require picking k upfront; the color threshold is tuned
   to roughly match the chosen k, for a consistent story across tabs.
9. **Confidence ellipses** (App only, optional): 80/85/90% confidence
   ellipses per group, on the PCA scatter plots, assuming a bivariate
   normal distribution -- a quick visual read on how tight or overlapping
   the groups/clusters are, beyond the raw scatter of points.

## Data

`Spectra/` (not tracked in this repository — see `.gitignore`) should
contain one `.dpt` file per spectrum: two comma-separated columns
(wavenumber, intensity), no header, same wavenumber grid across every
file in the folder.

## Included third-party/personal code

- `backcor.m` — V. Mazet et al., "Background removal from spectra by
  designing and minimising a non-quadratic cost function", Chemometrics
  and Intelligent Laboratory Systems, 76(2), 121-133 (2005).
  Implementation: V. Mazet (vincent.mazet@unistra.fr), public domain.
- `airPLS.m` — Zhang et al., "Baseline correction using adaptive
  iteratively reweighted penalized least squares", Analyst 135(5),
  1138-1146 (2010). Implementation: Zhimin Zhang, Central South
  University (2011), public domain (MATLAB File Exchange).
- `snip.m` — C. G. Ryan et al., "SNIP, a statistics-sensitive background
  treatment for the quantitative analysis of PIXE spectra in geoscience
  applications", Nuclear Instruments and Methods in Physics Research B,
  34, 396-402 (1988). Implementation: original, from the published
  algorithm description.
- `apls.m` — P. J. Cadusch et al., "Improved methods for fluorescence
  background subtraction from Raman spectra", J. Raman Spectrosc. 44(11),
  1587-1595 (2013). Implementation: original, from the published
  algorithm description.
- `readdpt.m` — personal library (`myfileutil/`), reads the plain
  two-column `.dpt` spectrum export format.

---

Developed with the assistance of an AI coding tool (Claude, Anthropic),
under the author's supervision and review.
