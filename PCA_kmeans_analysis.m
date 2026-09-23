%PCA_KMEANS_ANALYSIS  Unsupervised PCA + k-means clustering of the Raman
%   spectra in Spectra/.
%
%   Pipeline: load -> baseline/smooth/normalize -> PCA -> choose k
%   (elbow + silhouette over the PCA scores) -> k-means -> compare
%   clusters against the group each spectrum's filename suggests (only
%   used here for validation/plotting, never fed into the clustering
%   itself, which is unsupervised).
%
%   Figures and numeric results are saved under Results/.

spectraDir = fullfile(fileparts(mfilename('fullpath')), 'Spectra');
resultsDir = fullfile(fileparts(mfilename('fullpath')), 'Results');
if ~isfolder(resultsDir)
    mkdir(resultsDir);
end

%% 1. Load
[X, wavenumbers, labels, filenames] = loadRamanSpectra(spectraDir);
fprintf('Loaded %d spectra (%d points each, %.0f-%.0f cm^{-1}).\n', ...
    size(X,1), size(X,2), min(wavenumbers), max(wavenumbers));
groupNames = unique(labels);
fprintf('Filename-derived groups: %s\n', strjoin(groupNames, ', '));

%% 2. Preprocess (baseline removal, light smoothing, unit-area normalization)
[Xproc, baselines] = preprocessSpectra(X, wavenumbers);

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

%% 7. Visualize: PC1 vs PC2, colored by filename group and by cluster
fig = figure('Position', [100 100 1000 450]);
tiledlayout(fig, 1, 2);
nexttile;
gscatter(score(:,1), score(:,2), labels);
xlabel(sprintf('PC1 (%.1f%%)', explained(1))); ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title('Colored by filename-derived group');
nexttile;
gscatter(score(:,1), score(:,2), clusterIdx);
xlabel(sprintf('PC1 (%.1f%%)', explained(1))); ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title(sprintf('Colored by k-means cluster (k=%d)', k));
exportgraphics(fig, fullfile(resultsDir, 'pca_scatter.png'));

%% 8. Loadings (which wavenumbers drive PC1/PC2) and per-cluster mean spectra
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
        'DisplayName', sprintf('Cluster %d (n=%d)', c, sum(clusterIdx==c)));
end
xlabel('Raman shift (cm^{-1})'); ylabel('Normalized intensity (curves offset for clarity)');
set(gca, 'YTick', []);
legend('Location', 'best');
title('Mean spectrum per cluster');
exportgraphics(fig, fullfile(resultsDir, 'loadings_and_clusters.png'));

%% 9. Save numeric results
resultsTable = table(filenames, labels, clusterIdx, score(:,1), score(:,2), score(:,3), ...
    'VariableNames', {'FileName','FilenameGroup','Cluster','PC1','PC2','PC3'});
writetable(resultsTable, fullfile(resultsDir, 'cluster_assignments.csv'));
save(fullfile(resultsDir, 'pca_kmeans_results.mat'), 'X', 'Xproc', 'wavenumbers', 'labels', ...
    'filenames', 'coeff', 'score', 'explained', 'clusterIdx', 'contingency', 'ari', ...
    'kRange', 'wcss', 'meanSil', 'kSilhouette');

fprintf('\nDone. Figures and results saved to %s\n', resultsDir);
