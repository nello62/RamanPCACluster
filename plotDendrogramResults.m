function ax = plotDendrogramResults(source, varargin)
%PLOTDENDROGRAMRESULTS  Re-plot the Ward-linkage hierarchical clustering
%   dendrogram from a saved RamanPCACluster results file, without
%   rerunning the analysis.
%
%   AX = PLOTDENDROGRAMRESULTS(SOURCE) loads SOURCE (a .mat file path, or
%   an already-loaded struct -- see LOADPCARESULTS) and draws the
%   Ward-linkage dendrogram of the PCA scores that were actually
%   clustered on (SCOREREDUCED -- the same reduced subspace k-means used,
%   not just PC1-PC2), with its color threshold tuned to roughly match
%   the saved k -- the same plot as PCAClusterApp's own "Dendrogram" tab.
%
%   AX = PLOTDENDROGRAMRESULTS(..., 'Name', Value, ...) options:
%     'K'       target number of clusters for the color threshold
%               (default: derived from the saved CLUSTERIDX). Purely
%               cosmetic -- a dendrogram doesn't need k to be built, only
%               to decide where its coloring switches branches.
%     'Labels'  {n x 1} cellstr of leaf labels (default: FILENAMES, with
%               the file extension stripped).
%     'Axes'    target axes (default: creates a new figure/axes).
%
%   A results file saved before SCOREREDUCED was recorded falls back to
%   the full SCORE matrix instead, with a warning -- the dendrogram then
%   reflects every retained PC rather than just the ones k-means actually
%   clustered on, and may differ slightly from the one originally shown.
%
%   Example:
%     plotDendrogramResults('run1.mat');
%     plotDendrogramResults('run1.mat', 'K', 5);
%
%   See also LOADPCARESULTS, PLOTPCARESULTS, PLOTSCREERESULTS, LINKAGE,
%   DENDROGRAM, PCACLUSTERAPP.
    p = inputParser;
    addParameter(p, 'K', [], @(v) isempty(v) || (isscalar(v) && isnumeric(v) && v >= 1));
    addParameter(p, 'Labels', {}, @(v) iscellstr(v) || isstring(v));
    addParameter(p, 'Axes', [], @(v) isempty(v) || (isscalar(v) && isgraphics(v, 'axes')));
    parse(p, varargin{:});
    opt = p.Results;

    S = loadPCAResults(source);

    if isempty(S.scoreReduced)
        warning('plotDendrogramResults:noScoreReduced', ...
            ['This file predates saving scoreReduced (the exact PCA subspace k-means ' ...
             'clustered on) -- falling back to the full score matrix. The dendrogram may ' ...
             'differ slightly from the one originally shown.']);
        clusterSpace = S.score;
    else
        clusterSpace = S.scoreReduced;
    end

    if isempty(opt.K)
        k = max(S.clusterIdx);
    else
        k = round(opt.K);
    end

    if isempty(opt.Labels)
        [~, leafLabels] = cellfun(@fileparts, S.filenames, 'UniformOutput', false);
    else
        leafLabels = cellstr(opt.Labels);
    end

    Z = linkage(clusterSpace, 'ward', 'euclidean');
    nZ = size(Z, 1);
    if k >= 2 && k <= nZ
        cutoff = mean(Z(max(nZ - k + 1, 1):min(nZ - k + 2, nZ), 3));
    else
        cutoff = 0.7 * max(Z(:, 3));
    end

    if isempty(opt.Axes)
        fig = figure;
        ax = axes(fig);
    else
        ax = opt.Axes;
    end
    cla(ax);

    % DENDROGRAM predates axes-handle support and always draws into a
    % regular figure/GCA, silently ignoring any axes handle passed to it
    % (same limitation documented and worked around in PCACLUSTERAPP.M)
    % -- draw into a throwaway invisible figure instead, then copy the
    % resulting lines and tick setup into AX and discard the temporary
    % figure.
    tmpFig = figure('Visible', 'off');
    dendrogram(Z, 0, 'ColorThreshold', cutoff, 'Labels', leafLabels);
    tmpAx = gca;
    copyobj(tmpAx.Children, ax);
    ax.XTick = tmpAx.XTick;
    ax.XTickLabel = tmpAx.XTickLabel;
    ax.XLim = tmpAx.XLim;
    ax.YLim = tmpAx.YLim;
    close(tmpFig);

    ax.XTickLabelRotation = 90;
    xlabel(ax, 'Spectrum');
    ylabel(ax, 'Ward linkage distance');
    title(ax, sprintf('Hierarchical clustering dendrogram (color threshold tuned for k=%d)', k));
end
