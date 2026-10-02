# RamanPCACluster

Unsupervised PCA + k-means clustering of Raman spectra (`.dpt`/`.csv`/
`.spc`/`.wdf` format), as a MATLAB script and an interactive App.

## Requirements

- MATLAB with **Statistics and Machine Learning Toolbox** (`pca`, `kmeans`,
  `silhouette`, `fitgmdist`, `fitcdiscr`)
- **Signal Processing Toolbox** (`sgolayfilt`)
- No external dependencies beyond MATLAB itself: `backcor.m`, `airPLS.m`,
  `snip.m`, `apls.m`, `readdpt.m`, and the `.spc`/`.wdf` readers (see
  "Supported spectrum file formats" below) are included in this folder.

## Files

- `PCAClusterApp.m` — interactive app: pick a folder of spectra (any mix
  of `.dpt`/`.csv`/`.spc`/`.wdf`, see "Supported spectrum file formats"
  below), set the spectral range, preprocessing, and number of clusters, run the
  analysis, inspect the results across tabs (including a hierarchical
  clustering dendrogram, an LDA tab, a GMM soft-clustering tab, an
  Identification tab, and optional 80/85/90% confidence ellipses on the
  PCA scatter plots), and export them. "New session (clear all)" resets
  every setting/result back to the app's startup state to begin an
  unrelated analysis. A folder of known-material reference spectra (e.g.
  SLoPP/) can optionally be overlaid on the PCA scatter plots, one marker
  shape per material class -- purely for visual comparison, never fed
  into the analysis (see `loadReferenceSpectra.m` below) -- and used to
  suggest a material identity for each k-means cluster (see
  `identifyClustersByReference.m` below).
- `PCA_kmeans_analysis.m` — the same pipeline as a plain script (edit the
  top of the file to point at a different `Spectra/` folder), useful for
  batch/reproducible runs outside the GUI.
- `loadRamanSpectra.m` — loads every `.dpt`/`.csv`/`.spc`/`.wdf` file in a
  folder onto a common wavenumber grid (same readers as RamanFitApp, see
  "Supported spectrum file formats" below); also derives a group label
  per spectrum from its filename (everything before the underscore that
  introduces the sample identifier, e.g. `CT_P_11b.0.dpt` → `CT_P`,
  `Paradiso_31-nonpellet.0.dpt` → `Paradiso` — tolerant of
  irregularly-formatted sample identifiers), used only for
  validating/plotting the clustering, never for the clustering itself.
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
  `.txt`/`.dpt`/`.csv`/`.spc`/`.wdf` files, own native wavenumber grid per
  file -- unlike the main dataset, a reference library is not required to
  share one common grid, nor to all use the same format), and derives
  each one's material class from its file name (e.g. "Acrylic 1. Red
  Fiber (Polyacrylonitrile).txt" -> "Acrylic").
- `projectReferenceSpectra.m` — interpolates each reference spectrum onto
  the current analysis' own wavenumber grid, preprocesses it with the
  exact same options as the main dataset, and projects it into the
  already-fitted PCA space (reusing that fit's loadings/centering,
  never refitting PCA on the references). References whose own measured
  range does not fully cover the analysis range are skipped and reported,
  not silently dropped or extrapolated.
- `plotReferenceOverlay.m` — draws already-projected reference points on
  a PCA scatter (2D) or PCA scatter3 (3D, based on how many columns of
  score it's given) plot, one marker shape (cycling through a fixed set)
  and HSV color per material class, with a legend entry each.
- `identifyClustersByReference.m` — indicative material identification:
  for each k-means cluster, the reference material class whose spectra
  (averaged into that class's own centroid) sit closest in PCA space,
  plus the runner-up class/distance to gauge how confident or ambiguous
  the match is. A descriptive nearest-centroid heuristic, not a fitted or
  cross-validated classifier -- see "Cluster identification" in the
  Pipeline section below.
- `loadPCAResults.m`, `plotPCAResults.m`, `plotScreeResults.m`,
  `plotDendrogramResults.m`, `plotGMMResults.m`, `plotLDAResults.m`,
  `drawGaussianEllipse.m`, `drawGaussianEllipsoid.m` — re-plot a previously
  saved results file without rerunning the analysis; see "Re-plotting
  saved results" below.
- `backcor.m`, `airPLS.m`, `snip.m`, `apls.m`, `readdpt.m`, `readspc.m`,
  `readwdf.m` (plus the `GSSpcRead`/`WdfReader` third-party libraries they
  wrap -- see "Supported spectrum file formats" below) — third-party/
  personal helper functions, copied in for self-containment (see below).

## Pipeline

1. **Load**: every `.dpt`/`.csv`/`.spc`/`.wdf` spectrum in a folder,
   stacked into one matrix (all files must share the same wavenumber
   grid, regardless of format).
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
7.5. **Cluster identification** (indicative, only if a reference library
   is loaded): for each k-means cluster, the reference material class
   whose spectra sit closest by centroid distance in the same PCA
   subspace k-means clustered on, plus the runner-up class/distance --
   e.g. "Cluster 3: Polypropylene (d=0.0030); runner-up: Polyamide
   (d=0.0068)". A nearest-centroid heuristic, not a fitted/cross-validated
   classifier: it always returns *a* closest class, even from a small or
   incomplete library that may not actually contain the true material --
   read the distance (a small best-match distance with a much larger
   runner-up is a confident match; two nearly-tied distances is not) and
   the PCA scatter itself, not the label alone.
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

## Re-plotting saved results

`plotPCAResults.m`, `plotScreeResults.m`, `plotDendrogramResults.m`,
`plotGMMResults.m`, and `plotLDAResults.m` re-plot figures from a results
file saved by either app or script (`<base>.mat` / `pca_kmeans_results.mat`),
without rerunning any of the analysis -- a growing family of `plot*`
functions meant to work directly off a saved file.

```matlab
plotPCAResults('run1.mat');                          % PC1 vs PC2, k-means clusters
plotPCAResults('run1.mat', 'ColorBy', 'group');       % colored by filename-derived group instead
plotPCAResults('run1.mat', 'ColorBy', 'gmm');         % colored by GMM component (if the file has one)
plotPCAResults('run1.mat', 'Ellipses', false, 'Reference', false);  % bare scatter only
plotPCAResults('run1.mat', 'Dims', [1 2 3]);          % 3D PC1-PC2-PC3, with confidence ellipsoids

plotScreeResults('run1.mat');                         % explained + cumulative explained variance
plotScreeResults('run1.mat', 'NumComponents', 15);    % show more/fewer leading PCs than the default 10

plotDendrogramResults('run1.mat');                    % Ward-linkage dendrogram
plotDendrogramResults('run1.mat', 'K', 6);            % tune the color threshold for a different k

plotGMMResults('run1.mat');                           % GMM scatter (model ellipses) + BIC scan
plotGMMResults('run1.mat', 'Dims', [1 2 3]);          % 3D scatter with confidence ellipsoids

plotLDAResults('run1.mat');                           % LD1 vs LD2 + cross-validated confusion matrix
```

`plotPCAResults` draws the experimental data, 80/85/90% confidence
ellipses (2D) or ellipsoids (3D) per group/cluster, and the saved
reference-library overlay if the file has one -- the same information as
the app's own "Clusters" tab. `plotScreeResults` draws each PC's own
explained variance (bars) alongside the cumulative explained variance
(line) -- the same information as the app's own "PCA" tab.
`plotDendrogramResults` draws the Ward-linkage hierarchical clustering
tree on `scoreReduced` (the exact PCA subspace k-means clustered on, not
just PC1-PC2 -- saved alongside `score` since this function was added;
an older file without it falls back to the full `score` matrix with a
warning). `plotGMMResults` draws the GMM soft-clustering scatter (each
component's confidence ellipse/ellipsoid from the *fitted model's own*
mean/covariance, not a post-hoc empirical one) alongside the saved BIC
model-selection scan -- the same information as the app's own "GMM" tab;
errors if the file has no GMM result. `plotLDAResults` draws the LD1-vs-LD2
canonical discriminant scatter (colored by filename-derived group) and the
cross-validated confusion matrix -- the same information as the app's own
"LDA" tab; errors if the file has no LDA result (predates the feature, or
fewer than two filename-derived groups were present at analysis time). All
five reconstruct their figure from the `.mat` file alone. `loadPCAResults.m`
(loads/validates the file,
filling in any field an older save is missing) and
`drawGaussianEllipse.m`/`drawGaussianEllipsoid.m` (the 2D/3D
confidence-region primitives, shared with `PCAClusterApp.m` itself) are
shared building blocks, usable directly for a custom plot.

## Supported spectrum file formats

`loadRamanSpectra.m` scans a folder for `*.dpt`, `*.csv`, `*.spc`, and
`*.wdf` files (any mix of the four in the same folder is fine, as long as
every file shares the same wavenumber grid) and dispatches each one to
the matching reader -- the same readers RamanFitApp.m uses, so a folder
of spectra exported for that app can be reused here unchanged:

- **`.dpt`** — comma-separated wavenumber,intensity pairs, no header
  (`readdpt.m`; Hamamatsu/OPUS-style instrument export).
- **`.csv`** — exactly two numeric columns (wavenumber, intensity), no
  header row (plain `readmatrix`). Rejected, with a clear error naming
  the file, if it has a header row or more/fewer than two columns --
  this also catches a results `.csv` this app itself saved into the same
  folder (e.g. `<base>_cluster_assignments.csv`), rather than silently
  misreading it as a spectrum.
- **`.spc`** — Galactic/GRAMS binary spectrum format (`readspc.m`,
  wrapping `GSSpcRead.m`).
- **`.wdf`** — Renishaw WiRE binary spectrum format (`readwdf.m`,
  wrapping `WdfReader.m`).

Every file is expected to hold exactly one spectrum. Unlike RamanFitApp,
a "wide" multi-spectrum file (several intensity columns sharing one
wavenumber axis) is **not** supported here, since this loader assumes one
file = one row of the data matrix = one filename-derived group label.

The reference library (`loadReferenceSpectra.m`, "Sidebar: Reference
library") is a separate loader with its own, more permissive rules (no
shared-grid requirement, skips unreadable files with a warning instead of
erroring) -- it reads the same five formats (`.txt`/`.dpt`/`.csv`/`.spc`/
`.wdf`, any mix), described under "Data" below.

## Data

`Spectra/` (not tracked in this repository — see `.gitignore`) should
contain one spectrum file per file (`.dpt`/`.csv`/`.spc`/`.wdf`, see
"Supported spectrum file formats" above), all sharing the same
wavenumber grid across the folder.

`SLoPP/` (also not tracked — see `.gitignore`), when present, is used as
the default reference library: the Southern California Coastal Water
Research Project's plastics Raman reference spectra (SLoPP/SLoPP-E), one
`.txt` file per reference spectrum (two whitespace-separated columns,
wavenumber/intensity, no header, own native grid per file — unlike
`Spectra/`, these do not need to share a common grid). Any other folder
of reference spectra works too (`.txt`/`.dpt`/`.csv`/`.spc`/`.wdf`, any
mix -- picked via "Select reference folder..." in the app, or
`referenceDir` in the script); a non-spectrum file that happens to sit in
the folder (e.g. SLoPP's own `lib_names*.txt` index) is skipped with a
warning rather than breaking the load.

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
- `readspc.m` (wraps `GSSpcRead.m`, `GSSpcReadStructure.m`,
  `GetSPCAxisTypes.m`, `GetTechniques.m`, `GSToolsAbout.m`,
  `LocateItem.m`, `trimstr.m`, `trimleft.m`, `trimright.m`) — reads the
  Galactic Industries/Thermo GRAMS `.spc` binary spectrum format.
  Implementation: GSTools, by Kris De Gussem, Raman Spectroscopy Research
  Group, Department of Analytical Chemistry, Ghent University
  (2004-2009), dual-licensed GPLv3/BSD (full license text in each file's
  own header). Copied from `prog/GSTools/` (same subset RamanFitApp.m
  uses); `readspc.m` itself is original, a thin wrapper around
  `GSSpcRead`.
- `readwdf.m` (wraps `WdfReader.m`, `WdfError.m`, `WdfBlockID.m`,
  `WiREDataType.m`, `WiREDataUnit.m`, `WiREFocusMode.m`, `WiREKeys.m`,
  `WiREMeasurementType.m`, `WiREScanBasicType.m`, `WiREScanType.m`) —
  reads the Renishaw WiRE `.wdf` binary spectrum format, via Renishaw's
  own official MATLAB access package ("Renishaw WiRE WDF access package
  for MATLAB", MathWorks File Exchange / `Renishaw/wdf-matlab` on
  GitHub), copyright Renishaw plc, dual-licensed Apache-2.0/BSD-3-Clause.
  This subset covers read-only, single-spectrum access only (no write
  support); `readwdf.m` itself is original, a thin wrapper around
  `WdfReader`.

---

Developed with the assistance of an AI coding tool (Claude, Anthropic),
under the author's supervision and review.
