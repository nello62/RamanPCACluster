%PCA_KMEANS_ANALYSIS  Unsupervised PCA + k-means clustering of the Raman
%   spectra in Spectra/.
%
%   Pipeline: load -> baseline/smooth/normalize -> PCA -> choose k
%   (elbow + silhouette over the PCA scores) -> k-means -> compare
%   clusters against the group each spectrum's filename suggests (only
%   used here for validation/plotting, never fed into the clustering
%   itself, which is unsupervised) -> LDA (supervised check against the
%   filename-derived groups) -> Gaussian Mixture Model (soft-clustering
%   alternative to k-means, same k, via RUNGMM). A reference library
%   (REFERENCEDIR, e.g. the SLoPP/SLoPP-E plastics spectra) is optionally
%   projected into the same PCA space and overlaid on the scatter plots,
%   purely for visual comparison -- never fed into PCA/k-means/GMM/LDA.
%
%   Figures and numeric results are saved under Results/.

spectraDir = fullfile(fileparts(mfilename('fullpath')), 'Spectra');
resultsDir = fullfile(fileparts(mfilename('fullpath')), 'Results');
if ~isfolder(resultsDir)
    mkdir(resultsDir);
end

% Optional reference library (e.g. SLoPP/SLoPP-E known-material spectra),
% overlaid on the PCA scatter plots for visual comparison only -- set to
% '' to disable. Never enters PCA/k-means/GMM/LDA itself (see step 8.5).
referenceDir = fullfile(fileparts(mfilename('fullpath')), 'SLoPP');

%% 1. Load
[X, wavenumbers, labels, filenames] = loadRamanSpectra(spectraDir);
fprintf('Loaded %d spectra (%d points each, %.0f-%.0f cm^{-1}).\n', ...
    size(X,1), size(X,2), min(wavenumbers), max(wavenumbers));
groupNames = unique(labels);
fprintf('Filename-derived groups: %s\n', strjoin(groupNames, ', '));

%% 2. Preprocess (baseline removal, light smoothing, unit-area normalization)
% OPT is passed explicitly (rather than relying on PREPROCESSSPECTRA's own
% defaults implicitly, twice) so step 8.5's reference-spectra projection
% is guaranteed to use the exact same preprocessing as the main dataset.
opt = struct();
[Xproc, baselines] = preprocessSpectra(X, wavenumbers, opt);

fig = figure('Position', [100 100 900 400]);
tiledlayout(fig, 1, 2);
nexttile;
plot(wavenumbers, X(1,:), 'Color', [0.6 0.6 0.6]); hold on;
plot(wavenumbers, baselines(1,:), 'r--');
title('Example: raw spectrum + estimated baseline');
xlabel('Raman shift (cm^{-1})'); ylabel('Intensity');
legend({'Raw','Baseline'}, 'Location', 'best');
nexttile;
plot(wavenumbers, Xproc(1,:), 'b');
title('Example: after baseline/smooth/normalize'); xlabel('Raman shift (cm^{-1})'); ylabel('Normalized intensity');
exportgraphics(fig, fullfile(resultsDir, 'preprocessing_example.png'));

%% 3. PCA
[coeff, score, ~, ~, explained] = pca(Xproc);
nPCsFor95 = find(cumsum(explained) >= 95, 1, 'first');
fprintf('%d PCs needed to explain >=95%% of variance (PC1=%.1f%%, PC2=%.1f%%, PC3=%.1f%%).\n', ...
    nPCsFor95, explained(1), explained(2), explained(3));

fig = figure('Position', [100 100 700 400]);
yyaxis left
bar(explained(1:min(10,end)));
ylabel('Explained variance (%)');
yyaxis right
plot(cumsum(explained(1:min(10,end))), '-o');
ylabel('Cumulative explained variance (%)');
xlabel('Principal component');
title('Scree plot');
exportgraphics(fig, fullfile(resultsDir, 'scree_plot.png'));

% Cluster on enough PCs to explain most of the variance while dropping
% the noisiest high-order ones -- capped at 10 so a long tail of tiny
% components can't dilute the distances k-means uses.
nPCsForClustering = min(max(nPCsFor95, 2), 10);
scoreReduced = score(:, 1:nPCsForClustering);

%% 3.5. Reference library (optional): project known-material spectra into
% this PCA space for visual comparison. LOADREFERENCESPECTRA/
% PROJECTREFERENCESPECTRA never feed them into PCA/k-means/GMM/LDA --
% references skipped for not covering [min(wavenumbers) max(wavenumbers)]
% are reported, not silently dropped.
haveReferences = ~isempty(referenceDir) && isfolder(referenceDir);
if haveReferences
    [refWN, refIntensity, refClass, refNames] = loadReferenceSpectra(referenceDir);
    [refScore, refClassUsed, refNamesUsed, refSkipped] = projectReferenceSpectra( ...
        refWN, refIntensity, refClass, refNames, wavenumbers, opt, mean(Xproc, 1), coeff);
    fprintf('References: %d/%d used (%d classes), %d skipped (range not covered by %s-%s cm^{-1}).\n', ...
        numel(refNamesUsed), numel(refWN), numel(unique(refClassUsed)), numel(refSkipped), ...
        num2str(min(wavenumbers)), num2str(max(wavenumbers)));
else
    refScore = []; refClassUsed = {}; refNamesUsed = {}; refSkipped = {};
end

%% 4. Choose k: elbow (within-cluster sum of squares) + silhouette
kRange = 2:8;
wcss = zeros(size(kRange));
meanSil = zeros(size(kRange));
rng(1);  % reproducible k-means across runs of this script
for ik = 1:numel(kRange)
    k = kRange(ik);
    [idxK, ~, sumdK] = kmeans(scoreReduced, k, 'Replicates', 20);
    wcss(ik) = sum(sumdK);
    meanSil(ik) = mean(silhouette(scoreReduced, idxK));
end

fig = figure('Position', [100 100 900 400]);
tiledlayout(fig, 1, 2);
nexttile;
plot(kRange, wcss, '-o'); xlabel('k'); ylabel('Within-cluster sum of squares');
title('Elbow plot');
nexttile;
plot(kRange, meanSil, '-o'); xlabel('k'); ylabel('Mean silhouette');
title('Silhouette analysis');
exportgraphics(fig, fullfile(resultsDir, 'choosing_k.png'));

[~, bestKIdx] = max(meanSil);
kSilhouette = kRange(bestKIdx);
fprintf('Silhouette-suggested k = %d (mean silhouette = %.3f); %d filename-derived groups.\n', ...
    kSilhouette, meanSil(bestKIdx), numel(groupNames));

%% 5. Final k-means (k = number of filename-derived groups, for direct comparison)
k = numel(groupNames);
[clusterIdx, ~, ~] = kmeans(scoreReduced, k, 'Replicates', 50);

%% 6. Compare clusters against filename-derived groups
[groupIdx, ~] = grp2idx(labels);
contingency = crosstab(groupIdx, clusterIdx);
fprintf('\nContingency table (rows = filename-derived group, columns = k-means cluster):\n');
disp(array2table(contingency, 'VariableNames', strcat('Cluster', string(1:k)), 'RowNames', groupNames));

ari = adjustedRandIndex(groupIdx, clusterIdx);
fprintf('Adjusted Rand Index (agreement with filename-derived groups): %.3f\n', ari);

%% 6.5. Cluster identification (indicative): nearest reference material
% class per cluster, by centroid distance in the same PCA subspace
% k-means itself clustered on. A purely descriptive heuristic (nearest-
% centroid matching, no fitted/validated classifier) -- always read the
% distance alongside the label, not the label alone.
identBestClass = {}; identBestDist = []; identSecondClass = {}; identSecondDist = [];
if ~isempty(refScore)
    nDimsIdent = size(scoreReduced, 2);
    clusterCentroids = zeros(k, nDimsIdent);
    for c = 1:k
        clusterCentroids(c, :) = mean(scoreReduced(clusterIdx == c, :), 1);
    end
    [identBestClass, identBestDist, identSecondClass, identSecondDist] = ...
        identifyClustersByReference(clusterCentroids, refScore(:, 1:nDimsIdent), refClassUsed);
    fprintf('\nCluster identification (indicative, nearest reference class):\n');
    for c = 1:k
        fprintf('  Cluster %d (n=%d): %s (d=%.4f); runner-up: %s (d=%.4f)\n', ...
            c, sum(clusterIdx == c), identBestClass{c}, identBestDist(c), ...
            identSecondClass{c}, identSecondDist(c));
    end
    identTable = table((1:k)', arrayfun(@(c) sum(clusterIdx == c), (1:k)'), ...
        identBestClass(:), identBestDist(:), identSecondClass(:), identSecondDist(:), ...
        'VariableNames', {'Cluster','N','BestMatchClass','BestMatchDistance','RunnerUpClass','RunnerUpDistance'});
    writetable(identTable, fullfile(resultsDir, 'cluster_identification.csv'));
end

%% 7. LDA (supervised): how separable are the filename-derived groups
% themselves, independent of whatever k-means finds? Run on the same
% PCA-reduced scores used for clustering (LDA's within-class scatter
% matrix is singular on raw spectra, which have far more features than
% samples).
if numel(groupNames) > 1
    [ldaScores, explainedLDA, cvAccuracy, ldaConfMat, ldaClassNames] = runLDA(scoreReduced, labels);
    fprintf('LDA cross-validated classification accuracy (vs. filename-derived groups): %.1f%%\n', 100*cvAccuracy);

    fig = figure('Position', [100 100 1000 450]);
    tiledlayout(fig, 1, 2);
    nexttile;
    ldaColors = lines(numel(ldaClassNames));
    if size(ldaScores, 2) >= 2
        gscatter(ldaScores(:,1), ldaScores(:,2), labels, ldaColors);
        xlabel(sprintf('LD1 (%.1f%%)', explainedLDA(1))); ylabel(sprintf('LD2 (%.1f%%)', explainedLDA(2)));
    else
        gscatter((1:size(ldaScores,1))', ldaScores(:,1), labels, ldaColors);
        xlabel('Sample index'); ylabel(sprintf('LD1 (%.1f%%)', explainedLDA(1)));
    end
    title('LD1 vs LD2 (colored by filename-derived group)');
    nexttile;
    imagesc(ldaConfMat);
    colorbar;
    axis square;
    set(gca, 'XTick', 1:numel(ldaClassNames), 'XTickLabel', ldaClassNames, ...
        'YTick', 1:numel(ldaClassNames), 'YTickLabel', ldaClassNames, 'YDir', 'reverse');
    xlabel('Predicted (cross-validated)'); ylabel('Actual (filename-derived group)');
    for r = 1:size(ldaConfMat,1)
        for c2 = 1:size(ldaConfMat,2)
            text(c2, r, num2str(ldaConfMat(r,c2)), 'HorizontalAlignment', 'center', 'Color', 'w');
        end
    end
    title(sprintf('Confusion matrix (cross-validated accuracy = %.1f%%)', 100*cvAccuracy));
    exportgraphics(fig, fullfile(resultsDir, 'lda_results.png'));
else
    ldaScores = []; explainedLDA = []; cvAccuracy = []; ldaConfMat = []; ldaClassNames = {};
    fprintf('Only one filename-derived group found -- skipping LDA (needs at least two).\n');
end

%% 8. Gaussian Mixture Model: soft-clustering alternative to k-means
% Fitted at the same k as the final k-means above (for a direct,
% like-for-like comparison), plus a BIC scan over 1..8 components as a
% model-selection diagnostic in its own right.
[gmmClusterIdx, gmmPosterior, gmmBIC, gmmModel, gmmBICScan, gmmKScanUsed] = runGMM(scoreReduced, k, 1:8);
ariGMMvsKmeans = adjustedRandIndex(clusterIdx, gmmClusterIdx);
fprintf('GMM (k=%d) vs. k-means agreement (Adjusted Rand Index): %.3f\n', k, ariGMMvsKmeans);
if numel(groupNames) > 1
    ariGMMvsFilename = adjustedRandIndex(groupIdx, gmmClusterIdx);
    fprintf('GMM (k=%d) vs. filename-derived groups (Adjusted Rand Index): %.3f\n', k, ariGMMvsFilename);
else
    ariGMMvsFilename = NaN;
end

fig = figure('Position', [100 100 1000 450]);
tiledlayout(fig, 1, 2);
nexttile;
gmmColors = lines(size(gmmPosterior, 2));
gscatter(score(:,1), score(:,2), gmmClusterIdx, gmmColors);
xlabel(sprintf('PC1 (%.1f%%)', explained(1))); ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title(sprintf('GMM soft clustering (k=%d, mean max-posterior = %.2f)', ...
    size(gmmPosterior, 2), mean(max(gmmPosterior, [], 2))));
nexttile;
plot(gmmKScanUsed, gmmBICScan, '-o'); hold on;
[minBIC, ib] = min(gmmBICScan);
if ~isnan(minBIC)
    plot(gmmKScanUsed(ib), minBIC, 'r*', 'MarkerSize', 10);
end
yl = ylim; plot([k k], yl, 'k--');
xlabel('Number of components'); ylabel('BIC (lower is better)');
title(sprintf('Model selection: BIC vs. components (used k=%d, dashed line)', k));
exportgraphics(fig, fullfile(resultsDir, 'gmm_results.png'));

%% 9. Visualize: PC1 vs PC2, colored by filename group and by cluster
% Reference spectra (if any), when present, are overlaid on both panels
% with PLOTREFERENCEOVERLAY -- visual comparison only, never part of the
% clustering itself.
fig = figure('Position', [100 100 1000 450]);
tiledlayout(fig, 1, 2);
nexttile;
gscatter(score(:,1), score(:,2), labels);
xlabel(sprintf('PC1 (%.1f%%)', explained(1))); ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title('Colored by filename-derived group');
if ~isempty(refScore)
    hold on;
    plotReferenceOverlay(gca, refScore, refClassUsed);
    hold off;
end
nexttile;
clusterColors = lines(k);  % explicit, so the mean-spectra plot below can reuse the exact same colors
gscatter(score(:,1), score(:,2), clusterIdx, clusterColors);
xlabel(sprintf('PC1 (%.1f%%)', explained(1))); ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title(sprintf('Colored by k-means cluster (k=%d)', k));
if ~isempty(refScore)
    hold on;
    plotReferenceOverlay(gca, refScore, refClassUsed);
    hold off;
end
exportgraphics(fig, fullfile(resultsDir, 'pca_scatter.png'));

if ~isempty(refScore)
    refTable = table(refNamesUsed(:), refClassUsed(:), refScore(:,1), refScore(:,2), refScore(:,3), ...
        'VariableNames', {'Name','Class','PC1','PC2','PC3'});
    writetable(refTable, fullfile(resultsDir, 'reference_projections.csv'));
end

%% 10. Loadings (which wavenumbers drive PC1/PC2) and per-cluster mean spectra
% Both plotted as a vertical (waterfall-style) stack, each curve offset
% by a fixed step so overlapping peaks from different curves don't
% obscure each other; the Y axis ticks are hidden since the absolute
% offset position is then arbitrary, not a real intensity/loading value.
fig = figure('Position', [100 100 900 400]);
tiledlayout(fig, 1, 2);
nexttile;
l1 = coeff(:,1)';
l2 = coeff(:,2)';
loadingStep = 1.2 * max(range(l1), range(l2));
plot(wavenumbers, l1, 'DisplayName', 'PC1'); hold on;
plot(wavenumbers, l2 - loadingStep, 'DisplayName', 'PC2');
xlabel('Raman shift (cm^{-1})'); ylabel('Loading (curves offset for clarity)');
set(gca, 'YTick', []);
legend('Location', 'best');
title('PCA loadings');
nexttile;
meanSpectraAll = zeros(k, numel(wavenumbers));
for c = 1:k
    meanSpectraAll(c, :) = mean(Xproc(clusterIdx == c, :), 1);
end
clusterStep = 1.15 * max(range(meanSpectraAll, 2));
hold on;
for c = 1:k
    plot(wavenumbers, meanSpectraAll(c, :) + (c - 1) * clusterStep, ...
        'Color', clusterColors(c,:), 'DisplayName', sprintf('Cluster %d (n=%d)', c, sum(clusterIdx==c)));
end
xlabel('Raman shift (cm^{-1})'); ylabel('Normalized intensity (curves offset for clarity)');
set(gca, 'YTick', []);
legend('Location', 'best');
title('Mean spectrum per cluster');
exportgraphics(fig, fullfile(resultsDir, 'loadings_and_clusters.png'));

%% 11. Save numeric results
resultsTable = table(filenames, labels, clusterIdx, gmmClusterIdx, score(:,1), score(:,2), score(:,3), ...
    'VariableNames', {'FileName','FilenameGroup','KMeansCluster','GMMCluster','PC1','PC2','PC3'});
writetable(resultsTable, fullfile(resultsDir, 'cluster_assignments.csv'));
save(fullfile(resultsDir, 'pca_kmeans_results.mat'), 'X', 'Xproc', 'wavenumbers', 'labels', ...
    'filenames', 'coeff', 'score', 'explained', 'clusterIdx', 'contingency', 'ari', ...
    'kRange', 'wcss', 'meanSil', 'kSilhouette', ...
    'ldaScores', 'explainedLDA', 'cvAccuracy', 'ldaConfMat', 'ldaClassNames', ...
    'gmmClusterIdx', 'gmmPosterior', 'gmmModel', 'gmmBIC', 'gmmBICScan', 'gmmKScanUsed', ...
    'ariGMMvsKmeans', 'ariGMMvsFilename', ...
    'referenceDir', 'refScore', 'refClassUsed', 'refNamesUsed', 'refSkipped', ...
    'identBestClass', 'identBestDist', 'identSecondClass', 'identSecondDist');

fprintf('\nDone. Figures and results saved to %s\n', resultsDir);
