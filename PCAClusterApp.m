function PCAClusterApp()
%PCACLUSTERAPP  Interactive PCA + k-means clustering of a folder of Raman
%   spectra (.dpt format).
%
%   PCAClusterApp() opens the app; pick a folder of spectra, optionally
%   restrict the spectral range and preprocessing/clustering settings,
%   then press "Run analysis". Same underlying pipeline as
%   PCA_kmeans_analysis.m (LOADRAMANSPECTRA, PREPROCESSSPECTRA,
%   ADJUSTEDRANDINDEX), wrapped in a GUI so the folder, spectral range,
%   and number of clusters no longer require editing the script.
%
%   Requires the Statistics and Machine Learning Toolbox (PCA, KMEANS,
%   SILHOUETTE) and Signal Processing Toolbox (SGOLAYFILT, via
%   PREPROCESSSPECTRA).
%
%   See also LOADRAMANSPECTRA, PREPROCESSSPECTRA, ADJUSTEDRANDINDEX.

% -------------------------------------------------------------------------
% Session state (nested-function closures share these -- same single-file,
% no-classdef pattern used by RamanFitApp.m / G_gaussian_viewer.m).
% -------------------------------------------------------------------------
spectraDir = '';
X = []; wavenumbers = []; labels = {}; filenames = {};
Xproc = []; baselines = [];
coeff = []; score = []; explained = [];
clusterIdx = []; contingency = []; ari = []; k = [];
kRange = 2:8; wcss = []; meanSil = []; kSilhouette = [];
scoreReduced = [];
hasResults = false;

% -------------------------------------------------------------------------
fig = uifigure('Name', 'PCA Cluster App', 'Position', [100 100 1400 860]);

sidebarW = 340;
sidebar = uipanel(fig, 'Position', [0 0 sidebarW 860], 'BorderType', 'line');

statusLabel = uilabel(fig, 'Position', [sidebarW+10 6 1400-sidebarW-20 24], ...
    'Text', 'Select a folder of .dpt spectra to begin.', 'FontColor', [0.35 0.35 0.35]);

% ---- Sidebar: Data ------------------------------------------------------
uibutton(sidebar, 'push', 'Position', [10 820 sidebarW-20 30], ...
    'Text', 'Select spectra folder...', 'FontWeight', 'bold', ...
    'ButtonPushedFcn', @(s,e) onSelectFolder());
lblFolder = uilabel(sidebar, 'Position', [10 796 sidebarW-20 18], 'Text', 'Folder: -');
lblNSpectra = uilabel(sidebar, 'Position', [10 778 sidebarW-20 18], 'Text', 'Spectra: -');

% ---- Sidebar: Spectral range --------------------------------------------
uilabel(sidebar, 'Position', [10 744 sidebarW-20 18], 'Text', 'Spectral range (cm^{-1}):', 'FontWeight', 'bold');
uilabel(sidebar, 'Position', [10 718 30 18], 'Text', 'Min:');
rangeMinField = uieditfield(sidebar, 'numeric', 'Position', [42 716 110 22], 'Enable', 'off');
uilabel(sidebar, 'Position', [160 718 30 18], 'Text', 'Max:');
rangeMaxField = uieditfield(sidebar, 'numeric', 'Position', [192 716 110 22], 'Enable', 'off');
fullRangeBtn = uibutton(sidebar, 'push', 'Position', [10 684 sidebarW-20 26], ...
    'Text', 'Full range', 'Enable', 'off', 'ButtonPushedFcn', @(s,e) onFullRange());

% ---- Sidebar: Preprocessing ----------------------------------------------
uilabel(sidebar, 'Position', [10 648 sidebarW-20 18], 'Text', 'Preprocessing:', 'FontWeight', 'bold');
baselineCheck = uicheckbox(sidebar, 'Position', [10 622 sidebarW-20 22], ...
    'Text', 'Subtract baseline (airPLS)', 'Value', true, ...
    'ValueChangedFcn', @(s,e) onBaselineCheckChanged());
uilabel(sidebar, 'Position', [10 594 110 18], 'Text', 'Baseline lambda:');
baselineLambdaField = uieditfield(sidebar, 'numeric', 'Position', [140 592 sidebarW-160 22], 'Value', 1e6, 'Limits', [0 Inf]);
uilabel(sidebar, 'Position', [10 566 110 18], 'Text', 'Baseline order:');
baselineOrderField = uieditfield(sidebar, 'numeric', 'Position', [140 564 sidebarW-160 22], 'Value', 2, 'Limits', [1 4], 'RoundFractionalValues', 'on');
smoothCheck = uicheckbox(sidebar, 'Position', [10 536 sidebarW-20 22], 'Text', 'Savitzky-Golay smoothing', 'Value', true);
uilabel(sidebar, 'Position', [10 508 70 18], 'Text', 'Normalize:');
normalizeDD = uidropdown(sidebar, 'Position', [90 506 sidebarW-110 22], ...
    'Items', {'Area','Max','SNV','None'}, 'Value', 'Area');

% ---- Sidebar: Clustering --------------------------------------------------
uilabel(sidebar, 'Position', [10 472 sidebarW-20 18], 'Text', 'Clustering:', 'FontWeight', 'bold');
autoKCheck = uicheckbox(sidebar, 'Position', [10 446 sidebarW-20 22], ...
    'Text', 'Choose k automatically (silhouette, 2-8)', 'Value', false, ...
    'ValueChangedFcn', @(s,e) onAutoKChanged());
uilabel(sidebar, 'Position', [10 418 100 18], 'Text', 'Number of clusters k:');
kField = uieditfield(sidebar, 'numeric', 'Position', [140 416 sidebarW-160 22], 'Value', 3, 'Limits', [2 20], 'RoundFractionalValues', 'on');
ellipseCheck = uicheckbox(sidebar, 'Position', [10 386 sidebarW-20 22], ...
    'Text', 'Show confidence ellipses (80/85/90%)', 'Value', false, ...
    'ValueChangedFcn', @(s,e) onEllipseCheckChanged());

runBtn = uibutton(sidebar, 'push', 'Position', [10 344 sidebarW-20 34], ...
    'Text', 'Run analysis', 'FontWeight', 'bold', 'Enable', 'off', ...
    'ButtonPushedFcn', @(s,e) onRunAnalysis());
% Separate from STATUSLABEL (which reports the more detailed step-by-step
% progress at the bottom of the window): this one sits right under the
% button so it's impossible to miss, especially useful when re-running
% after just tweaking a parameter -- otherwise the plots still showing
% the PREVIOUS run's results could easily be mistaken for the new ones
% while the (possibly slow) computation is still in progress.
runningLabel = uilabel(sidebar, 'Position', [10 316 sidebarW-20 18], ...
    'Text', '', 'FontWeight', 'bold', 'FontColor', [0.85 0.35 0], 'HorizontalAlignment', 'center');

% ---- Sidebar: Export -------------------------------------------------------
uilabel(sidebar, 'Position', [10 274 sidebarW-20 18], 'Text', 'Export:', 'FontWeight', 'bold');
saveBtn = uibutton(sidebar, 'push', 'Position', [10 244 sidebarW-20 28], ...
    'Text', 'Save results...', 'Enable', 'off', 'ButtonPushedFcn', @(s,e) onSaveResults());

% ---- Main area: tabbed results ---------------------------------------------
tg = uitabgroup(fig, 'Position', [sidebarW+10 40 1400-sidebarW-20 810]);
tabPreprocess = uitab(tg, 'Title', 'Preprocessing');
tabPCA = uitab(tg, 'Title', 'PCA');
tabChooseK = uitab(tg, 'Title', 'Choose k');
tabClusters = uitab(tg, 'Title', 'Clusters');
tabDendro = uitab(tg, 'Title', 'Dendrogram');
tabLoadings = uitab(tg, 'Title', 'Loadings');

axPre1 = uiaxes(tabPreprocess, 'Position', [10 10 490 760]);
title(axPre1, 'Example: raw spectrum + estimated baseline');
axPre2 = uiaxes(tabPreprocess, 'Position', [510 10 490 760]);
title(axPre2, 'Example: after baseline/smooth/normalize');

axScree = uiaxes(tabPCA, 'Position', [10 10 990 760]);
title(axScree, 'Scree plot');

axElbow = uiaxes(tabChooseK, 'Position', [10 10 490 760]);
title(axElbow, 'Elbow plot');
axSil = uiaxes(tabChooseK, 'Position', [510 10 490 760]);
title(axSil, 'Silhouette analysis');

axClustersByGroup = uiaxes(tabClusters, 'Position', [10 10 490 760]);
title(axClustersByGroup, 'Colored by filename-derived group');
axClustersByKmeans = uiaxes(tabClusters, 'Position', [510 10 490 760]);
title(axClustersByKmeans, 'Colored by k-means cluster');

axDendro = uiaxes(tabDendro, 'Position', [10 10 990 760]);
title(axDendro, 'Hierarchical clustering dendrogram');

axLoadings = uiaxes(tabLoadings, 'Position', [10 10 490 760]);
title(axLoadings, 'PCA loadings');
axMeanSpectra = uiaxes(tabLoadings, 'Position', [510 10 490 760]);
title(axMeanSpectra, 'Mean spectrum per cluster');

% =========================================================================
%  Callbacks
% =========================================================================

    function clearAxesFully(ax)
    % CLA does not delete children whose HANDLEVISIBILITY is 'off' --
    % confirmed empirically here (and previously in RamanFitApp): the
    % confidence-ellipse lines in PLOTCLUSTERS are drawn with
    % HandleVisibility off precisely so they don't clutter the legend,
    % which means a plain CLA() before each re-plot left them stranded,
    % re-accumulating on every re-run with a new k. FINDALL, unlike CLA
    % or FINDOBJ, ignores HandleVisibility, so this actually removes
    % everything; the axes itself is filtered back out of its own
    % FINDALL results first so it doesn't get deleted along with them.
        kids = findall(ax);
        kids(kids == ax) = [];
        delete(kids);
    end

% -------------------------------------------------------------------------
    function onSelectFolder()
        d = uigetdir(pwd, 'Select a folder of .dpt Raman spectra');
        if isequal(d, 0)
            return
        end
        % Reuses RUNNINGLABEL (see ONRUNANALYSIS) rather than a separate
        % label: the two are never active at the same time, and both
        % exist for the same reason -- a prominent "something is
        % happening" cue, since a folder with many spectra can take a
        % moment to read.
        runningLabel.Text = 'Importing files...';
        drawnow;
        try
            [X, wavenumbers, labels, filenames] = loadRamanSpectra(d);
        catch ME
            runningLabel.Text = '';
            uialert(fig, ME.message, 'Load error');
            return
        end
        runningLabel.Text = '';
        spectraDir = d;
        [~, folderName] = fileparts(d);
        lblFolder.Text = sprintf('Folder: %s', folderName);
        lblNSpectra.Text = sprintf('Spectra: %d', size(X, 1));
        rangeMinField.Value = min(wavenumbers);
        rangeMaxField.Value = max(wavenumbers);
        rangeMinField.Enable = 'on';
        rangeMaxField.Enable = 'on';
        fullRangeBtn.Enable = 'on';
        runBtn.Enable = 'on';
        hasResults = false;
        saveBtn.Enable = 'off';
        statusLabel.Text = sprintf('Loaded %d spectra from %s. Set the range/options and press Run analysis.', ...
            size(X, 1), folderName);
    end

% -------------------------------------------------------------------------
    function onFullRange()
        if isempty(wavenumbers)
            return
        end
        rangeMinField.Value = min(wavenumbers);
        rangeMaxField.Value = max(wavenumbers);
    end

% -------------------------------------------------------------------------
    function onAutoKChanged()
        if autoKCheck.Value
            kField.Enable = 'off';
        else
            kField.Enable = 'on';
        end
    end

% -------------------------------------------------------------------------
    function onBaselineCheckChanged()
        if baselineCheck.Value
            baselineLambdaField.Enable = 'on';
            baselineOrderField.Enable = 'on';
        else
            baselineLambdaField.Enable = 'off';
            baselineOrderField.Enable = 'off';
        end
    end

% -------------------------------------------------------------------------
    function onEllipseCheckChanged()
    % Cheap enough to just redraw from the results already in memory --
    % no need to rerun PCA/k-means for a purely cosmetic toggle.
        if hasResults
            plotClusters(k);
        end
    end

% -------------------------------------------------------------------------
    function onRunAnalysis()
        if isempty(X)
            return
        end
        rangeMin = rangeMinField.Value;
        rangeMax = rangeMaxField.Value;
        if rangeMin >= rangeMax
            uialert(fig, 'Min must be less than Max.', 'Invalid range');
            return
        end
        mask = wavenumbers >= rangeMin & wavenumbers <= rangeMax;
        if nnz(mask) < 10
            uialert(fig, 'Selected range is too narrow (fewer than 10 points).', 'Invalid range');
            return
        end

        runBtn.Enable = 'off';
        saveBtn.Enable = 'off';
        runningLabel.Text = 'Running analysis...';
        statusLabel.Text = 'Running: preprocessing...';
        drawnow;

        wnRange = wavenumbers(mask);
        Xrange = X(:, mask);

        opt = struct('Baseline', logical(baselineCheck.Value), ...
            'BaselineLambda', baselineLambdaField.Value, ...
            'BaselineOrder', baselineOrderField.Value, ...
            'Smooth', logical(smoothCheck.Value), ...
            'Normalize', lower(normalizeDD.Value));
        try
            [Xproc, baselines] = preprocessSpectra(Xrange, wnRange, opt);
        catch ME
            uialert(fig, ME.message, 'Preprocessing error');
            runningLabel.Text = '';
            runBtn.Enable = 'on';
            if hasResults
                saveBtn.Enable = 'on';
            end
            return
        end
        wavenumbersUsed = wnRange;

        plotPreprocessingExample(Xrange, wavenumbersUsed);

        statusLabel.Text = 'Running: PCA...';
        drawnow;
        [coeff, score, ~, ~, explained] = pca(Xproc);
        plotScree();

        nPCsFor95 = find(cumsum(explained) >= 95, 1, 'first');
        nPCsForClustering = min(max(nPCsFor95, 2), 10);
        scoreReduced = score(:, 1:nPCsForClustering);

        statusLabel.Text = 'Running: choosing k...';
        drawnow;
        rng(1);
        wcss = zeros(size(kRange));
        meanSil = zeros(size(kRange));
        for ik = 1:numel(kRange)
            kk = kRange(ik);
            [idxK, ~, sumdK] = kmeans(scoreReduced, kk, 'Replicates', 20);
            wcss(ik) = sum(sumdK);
            meanSil(ik) = mean(silhouette(scoreReduced, idxK));
        end
        [~, bestKIdx] = max(meanSil);
        kSilhouette = kRange(bestKIdx);
        plotChooseK();

        if autoKCheck.Value
            k = kSilhouette;
            kField.Value = k;
        else
            k = kField.Value;
        end

        statusLabel.Text = sprintf('Running: k-means (k=%d)...', k);
        drawnow;
        clusterIdx = kmeans(scoreReduced, k, 'Replicates', 50);

        groupNames = unique(labels);
        if numel(groupNames) > 1
            [groupIdx, ~] = grp2idx(labels);
            contingency = crosstab(groupIdx, clusterIdx);
            ari = adjustedRandIndex(groupIdx, clusterIdx);
        else
            contingency = [];
            ari = NaN;
        end

        plotClusters(k);
        plotDendrogram(k);
        plotLoadingsAndClusters(k, wavenumbersUsed);

        hasResults = true;
        saveBtn.Enable = 'on';
        runBtn.Enable = 'on';
        runningLabel.Text = '';
        if isnan(ari)
            statusLabel.Text = sprintf('Done: k=%d clusters (silhouette-suggested k=%d).', k, kSilhouette);
        else
            statusLabel.Text = sprintf('Done: k=%d clusters (silhouette-suggested k=%d); Adjusted Rand Index vs. filename groups = %.3f.', ...
                k, kSilhouette, ari);
        end
    end

% -------------------------------------------------------------------------
    function plotPreprocessingExample(Xrange, wn)
        clearAxesFully(axPre1); clearAxesFully(axPre2);
        plot(axPre1, wn, Xrange(1,:), 'Color', [0.6 0.6 0.6]);
        if baselineCheck.Value
            hold(axPre1, 'on');
            plot(axPre1, wn, baselines(1,:), 'r--');
            hold(axPre1, 'off');
            legend(axPre1, {'Raw','Baseline'}, 'Location', 'best');
            title(axPre1, sprintf('Example (%s): raw spectrum + estimated baseline', filenames{1}));
        else
            legend(axPre1, {'Raw'}, 'Location', 'best');
            title(axPre1, sprintf('Example (%s): raw spectrum (baseline subtraction off)', filenames{1}));
        end
        xlabel(axPre1, 'Raman shift (cm^{-1})'); ylabel(axPre1, 'Intensity');

        plot(axPre2, wn, Xproc(1,:), 'b');
        xlabel(axPre2, 'Raman shift (cm^{-1})'); ylabel(axPre2, 'Normalized intensity');
        title(axPre2, 'Example: after baseline/smooth/normalize');
    end

% -------------------------------------------------------------------------
    function plotScree()
        clearAxesFully(axScree);
        nShow = min(10, numel(explained));
        yyaxis(axScree, 'left');
        bar(axScree, explained(1:nShow));
        ylabel(axScree, 'Explained variance (%)');
        yyaxis(axScree, 'right');
        plot(axScree, cumsum(explained(1:nShow)), '-o');
        ylabel(axScree, 'Cumulative explained variance (%)');
        xlabel(axScree, 'Principal component');
        title(axScree, 'Scree plot');
    end

% -------------------------------------------------------------------------
    function plotChooseK()
        clearAxesFully(axElbow); clearAxesFully(axSil);
        plot(axElbow, kRange, wcss, '-o');
        xlabel(axElbow, 'k'); ylabel(axElbow, 'Within-cluster sum of squares');
        title(axElbow, 'Elbow plot');
        plot(axSil, kRange, meanSil, '-o');
        xlabel(axSil, 'k'); ylabel(axSil, 'Mean silhouette');
        title(axSil, sprintf('Silhouette analysis (best k=%d)', kSilhouette));
    end

% -------------------------------------------------------------------------
    function plotClusters(k)
        clearAxesFully(axClustersByGroup); clearAxesFully(axClustersByKmeans);
        hG = gscatter(axClustersByGroup, score(:,1), score(:,2), labels);
        xlabel(axClustersByGroup, sprintf('PC1 (%.1f%%)', explained(1)));
        ylabel(axClustersByGroup, sprintf('PC2 (%.1f%%)', explained(2)));
        title(axClustersByGroup, 'Colored by filename-derived group');
        if ellipseCheck.Value
            groupNames = unique(labels);
            hold(axClustersByGroup, 'on');
            for gi = 1:numel(groupNames)
                gmask = strcmp(labels, groupNames{gi});
                plotConfidenceEllipse(axClustersByGroup, score(gmask,1), score(gmask,2), hG(gi).Color);
            end
            hold(axClustersByGroup, 'off');
        end

        hK = gscatter(axClustersByKmeans, score(:,1), score(:,2), clusterIdx);
        xlabel(axClustersByKmeans, sprintf('PC1 (%.1f%%)', explained(1)));
        ylabel(axClustersByKmeans, sprintf('PC2 (%.1f%%)', explained(2)));
        title(axClustersByKmeans, sprintf('Colored by k-means cluster (k=%d)', k));
        if ellipseCheck.Value
            hold(axClustersByKmeans, 'on');
            for c = 1:k
                cmask = clusterIdx == c;
                plotConfidenceEllipse(axClustersByKmeans, score(cmask,1), score(cmask,2), hK(c).Color);
            end
            hold(axClustersByKmeans, 'off');
        end
    end

% -------------------------------------------------------------------------
    function plotConfidenceEllipse(ax, x, y, color)
    % Draws the 80/85/90% confidence ellipses (dotted/dashed/solid) for a
    % single group of 2D points, assuming a bivariate normal distribution
    % -- a standard way to visualize how tight/overlapping clusters are
    % in PC1-PC2 space, beyond just the scatter of points itself.
    % HANDLEVISIBILITY off so these don't add clutter entries to the
    % existing per-group/per-cluster legend from GSCATTER.
        if numel(x) < 3
            return
        end
        mu = [mean(x), mean(y)];
        C = cov(x, y);
        [V, D] = eig(C);
        theta = linspace(0, 2*pi, 100);
        circle = [cos(theta); sin(theta)];
        confLevels = [0.80 0.85 0.90];
        styles = {':', '--', '-'};
        for i = 1:numel(confLevels)
            r = sqrt(chi2inv(confLevels(i), 2));
            pts = mu' + V * sqrt(D) * r * circle;
            plot(ax, pts(1,:), pts(2,:), styles{i}, 'Color', color, ...
                'LineWidth', 1.2, 'HandleVisibility', 'off');
        end
    end

% -------------------------------------------------------------------------
    function plotDendrogram(k)
    % Hierarchical clustering (Ward linkage) on the same PCA scores used
    % for k-means -- an alternative view of cluster structure that
    % doesn't require picking k upfront; COLORTHRESHOLD is tuned so the
    % dendrogram's own coloring lines up with the current k for a
    % consistent story across tabs, not because k literally means
    % anything to a dendrogram on its own.
        clearAxesFully(axDendro);
        Z = linkage(scoreReduced, 'ward', 'euclidean');
        nZ = size(Z, 1);
        if k >= 2 && k <= nZ
            cutoff = mean(Z(max(nZ-k+1,1):min(nZ-k+2,nZ), 3));
        else
            cutoff = 0.7 * max(Z(:,3));
        end
        leafLabels = regexprep(filenames, '\.dpt$', '');

        % DENDROGRAM predates UIAXES support and always draws into a
        % regular figure/GCA, silently ignoring any axes handle passed to
        % it (confirmed empirically: no error, just nothing drawn where
        % expected) -- draw into a throwaway invisible figure instead,
        % then copy the resulting lines and tick setup into AXDENDRO and
        % discard the temporary figure.
        tmpFig = figure('Visible', 'off');
        dendrogram(Z, 0, 'ColorThreshold', cutoff, 'Labels', leafLabels);
        tmpAx = gca;
        copyobj(tmpAx.Children, axDendro);
        axDendro.XTick = tmpAx.XTick;
        axDendro.XTickLabel = tmpAx.XTickLabel;
        axDendro.XLim = tmpAx.XLim;
        axDendro.YLim = tmpAx.YLim;
        close(tmpFig);

        axDendro.XTickLabelRotation = 90;
        xlabel(axDendro, 'Spectrum');
        ylabel(axDendro, 'Ward linkage distance');
        title(axDendro, sprintf('Hierarchical clustering dendrogram (color threshold tuned for k=%d)', k));
    end

% -------------------------------------------------------------------------
    function plotLoadingsAndClusters(k, wn)
    % Both PCA loadings and per-cluster mean spectra are plotted as a
    % vertical (waterfall-style) stack, each curve offset by a fixed
    % step so overlapping peaks from different curves don't obscure each
    % other -- standard practice for comparing multiple spectra/loadings
    % at a glance. The absolute Y position is then meaningless (only
    % arbitrary offsets), so the Y axis ticks are hidden rather than
    % showing numbers that don't mean anything on their own.
        clearAxesFully(axLoadings); clearAxesFully(axMeanSpectra);

        l1 = coeff(:,1)';
        l2 = coeff(:,2)';
        loadingStep = 1.2 * max(range(l1), range(l2));
        hold(axLoadings, 'on');
        plot(axLoadings, wn, l1, 'DisplayName', 'PC1');
        plot(axLoadings, wn, l2 - loadingStep, 'DisplayName', 'PC2');
        hold(axLoadings, 'off');
        legend(axLoadings, 'Location', 'best');
        xlabel(axLoadings, 'Raman shift (cm^{-1})');
        ylabel(axLoadings, 'Loading (curves offset for clarity)');
        axLoadings.YTick = [];
        title(axLoadings, 'PCA loadings');

        meanSpectraAll = zeros(k, numel(wn));
        for c = 1:k
            meanSpectraAll(c, :) = mean(Xproc(clusterIdx == c, :), 1);
        end
        clusterStep = 1.15 * max(range(meanSpectraAll, 2));
        hold(axMeanSpectra, 'on');
        for c = 1:k
            plot(axMeanSpectra, wn, meanSpectraAll(c, :) + (c - 1) * clusterStep, ...
                'DisplayName', sprintf('Cluster %d (n=%d)', c, sum(clusterIdx == c)));
        end
        hold(axMeanSpectra, 'off');
        legend(axMeanSpectra, 'Location', 'best');
        xlabel(axMeanSpectra, 'Raman shift (cm^{-1})');
        ylabel(axMeanSpectra, 'Normalized intensity (curves offset for clarity)');
        axMeanSpectra.YTick = [];
        title(axMeanSpectra, 'Mean spectrum per cluster');
    end

% -------------------------------------------------------------------------
    function onSaveResults()
        if ~hasResults
            return
        end
        d = uigetdir(spectraDir, 'Select a folder to save results into');
        if isequal(d, 0)
            return
        end
        exportgraphics(axPre1, fullfile(d, 'preprocessing_raw_and_baseline.png'));
        exportgraphics(axPre2, fullfile(d, 'preprocessing_processed.png'));
        exportgraphics(axScree, fullfile(d, 'scree_plot.png'));
        exportgraphics(axElbow, fullfile(d, 'elbow_plot.png'));
        exportgraphics(axSil, fullfile(d, 'silhouette_plot.png'));
        exportgraphics(axClustersByGroup, fullfile(d, 'clusters_by_filename_group.png'));
        exportgraphics(axClustersByKmeans, fullfile(d, 'clusters_by_kmeans.png'));
        exportgraphics(axDendro, fullfile(d, 'dendrogram.png'));
        exportgraphics(axLoadings, fullfile(d, 'pca_loadings.png'));
        exportgraphics(axMeanSpectra, fullfile(d, 'mean_spectrum_per_cluster.png'));

        resultsTable = table(filenames, labels, clusterIdx, score(:,1), score(:,2), score(:,3), ...
            'VariableNames', {'FileName','FilenameGroup','Cluster','PC1','PC2','PC3'});
        writetable(resultsTable, fullfile(d, 'cluster_assignments.csv'));
        save(fullfile(d, 'pca_kmeans_results.mat'), 'X', 'Xproc', 'wavenumbers', 'labels', ...
            'filenames', 'coeff', 'score', 'explained', 'clusterIdx', 'contingency', 'ari', ...
            'kRange', 'wcss', 'meanSil', 'kSilhouette');

        statusLabel.Text = sprintf('Results saved to %s.', d);
    end

end
