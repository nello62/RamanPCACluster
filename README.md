# RamanPCACluster

Unsupervised PCA + k-means clustering of Raman spectra (`.dpt` format),
as a MATLAB script and an interactive App.

## Requirements

- MATLAB with **Statistics and Machine Learning Toolbox** (`pca`, `kmeans`,
  `silhouette`, `fitgmdist`, `fitcdiscr`)
- **Signal Processing Toolbox** (`sgolayfilt`)
- No external dependencies beyond MATLAB itself: `backcor.m`, `airPLS.m`,
  `snip.m`, `apls.m`, and `readdpt.m` are included in this folder.

## Files

- `PCAClusterApp.m` — interactive app: pick a folder of `.dpt` spectra,
  set the spectral range, preprocessing, and number of clusters, run the
  analysis, inspect the results across tabs (including a hierarchical
  clustering dendrogram, an LDA tab, a GMM soft-clustering tab, and
  optional 80/85/90% confidence ellipses on the PCA scatter plots), and
  export them. A folder of known-material reference spectra (e.g. SLoPP/)
  can optionally be overlaid on the PCA scatter plots, one marker shape
  per material class -- purely for visual comparison, never fed into the
  analysis (see `loadReferenceSpectra.m` below).
- `PCA_kmeans_analysis.m` — the same pipeline as a plain script (edit the
  top of the file to point at a different `Spectra/` folder), useful for
  batch/reproducible runs outside the GUI.
- `loadRamanSpectra.m` — loads every `.dpt` file in a folder onto a common
  wavenumber grid; also derives a group label per spectrum from its
  filename (everything before the underscore that introduces the sample
  identifier, e.g. `CT_P_11b.0.dpt` → `CT_P`, `Paradiso_31-nonpellet.0.dpt`
  → `Paradiso` — tolerant of irregularly-formatted sample identifiers),
  used only for validating/plotting the clustering, never for the
  clustering itself.
- `preprocessSpectra.m` — baseline removal (choice of backcor / airPLS /
  SNIP / APLS, optional), Savitzky-Golay smoothing (optional), and
  normalization (area / max / SNV / none) -- the same four baseline
  methods RamanFitApp offers.
- `adjustedRandIndex.m` — Adjusted Rand Index between two partitions
  (e.g. k-means clusters vs. filename-derived groups), corrected for
  chance agreement.
- `runLDA.m` — supervised Linear Discriminant Analysis against the
  filename-derived groups: canonical discriminant scores (LD1, LD2, ...)
  for visualization, plus k-fold cross-validated classification accuracy
  and confusion matrix, as a check on how separable the groups actually
  are (independent of whatever k-means finds).
- `runGMM.m` — unsupervised Gaussian Mixture Model (Expectation-
  Maximization), fitted on the same PCA scores and at the same k as
  k-means, as a soft-clustering alternative: each spectrum gets a full
  posterior probability of belonging to every component instead of a
  single hard label, and each component is a full (not necessarily
  round) Gaussian in PCA-score space. Also scans 1-8 components and
  reports each fit's BIC, as a model-selection diagnostic independent of
  k-means' own elbow/silhouette choice of k.
- `loadReferenceSpectra.m` — loads a folder of reference spectra (flat
  `.txt` files, two whitespace-separated columns, own native wavenumber
  grid per file -- unlike the main dataset, a reference library is not
  required to share one common grid), and derives each one's material
  class from its file name (e.g. "Acrylic 1. Red Fiber
  (Polyacrylonitrile).txt" -> "Acrylic").
- `projectReferenceSpectra.m` — interpolates each reference spectrum onto
  the current analysis' own wavenumber grid, preprocesses it with the
  exact same options as the main dataset, and projects it into the
  already-fitted PCA space (reusing that fit's loadings/centering,
  never refitting PCA on the references). References whose own measured
  range does not fully cover the analysis range are skipped and reported,
  not silently dropped or extrapolated.
- `plotReferenceOverlay.m` — draws already-projected reference points on
  a PCA scatter plot, one marker shape (cycling through a fixed set) and
  HSV color per material class, with a legend entry each.
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
4.5. **Reference library** (optional): known-material reference spectra
   are interpolated onto the analysis grid, preprocessed identically to
   the main dataset, and projected into this same PCA space -- for visual
   comparison on the PCA scatter plots only, never fed into PCA/k-means/
   GMM/LDA. References whose own measured range doesn't cover the
   analysis range are skipped and reported.
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
10. **LDA** (supervised, optional check): using the filename-derived
    groups as known labels, canonical discriminant scores (LD1 vs. LD2)
    plus k-fold cross-validated classification accuracy and confusion
    matrix -- a way to ask "are these groups actually separable at all
    (by *any* method, not just unsupervised k-means)?", independent of
    the clustering itself. Needs at least two filename-derived groups.
11. **GMM** (unsupervised, soft-clustering alternative to k-means):
    fitted at the same k, on the same PCA scores, via Expectation-
    Maximization. Unlike k-means' hard, equal-size-favoring partition,
    every spectrum gets a posterior probability of belonging to each
    component, and each component's covariance is fitted freely rather
    than assumed spherical -- so it can capture elongated/correlated
    clusters k-means would split or merge incorrectly. Also reports the
    Adjusted Rand Index between the GMM and k-means partitions (how much
    the two methods agree), and a BIC scan over 1-8 components as an
    independent check on the chosen k.

## Data

`Spectra/` (not tracked in this repository — see `.gitignore`) should
contain one `.dpt` file per spectrum: two comma-separated columns
(wavenumber, intensity), no header, same wavenumber grid across every
file in the folder.

`SLoPP/` (also not tracked — see `.gitignore`), when present, is used as
the default reference library: the Southern California Coastal Water
Research Project's plastics Raman reference spectra (SLoPP/SLoPP-E), one
`.txt` file per reference spectrum (two whitespace-separated columns,
wavenumber/intensity, no header, own native grid per file — unlike
`Spectra/`, these do not need to share a common grid). Any other folder
of `.txt` reference spectra in the same format works too (picked via
"Select reference folder..." in the app, or `referenceDir` in the script).

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
