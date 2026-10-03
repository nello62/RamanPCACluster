function ax = plotIdentification3DResults(source, varargin)
%PLOTIDENTIFICATION3DRESULTS  Re-plot the PC1-PC2-PC3 cluster-centroid-to-
%   reference-class-centroid view from a saved RamanPCACluster results
%   file, without rerunning the analysis.
%
%   AX = PLOTIDENTIFICATION3DRESULTS(SOURCE) loads SOURCE (a .mat file
%   path, or an already-loaded struct -- see LOADPCARESULTS) and draws
%   into a new figure: a PC1-PC2-PC3 scatter of the samples (by k-means
%   cluster), the reference-library overlay, and a dashed line, per
%   cluster, from its own centroid (filled circle) to its matched
%   reference class's centroid (dark diamond), labeled with the class
%   name and distance in the legend -- the same information as
%   PCAClusterApp's own "Identification (3D)" tab, reconstructed from the
%   saved file alone.
%
%   AX = PLOTIDENTIFICATION3DRESULTS(..., 'Name', Value, ...) options:
%     'Ellipses'  true (default) | false -- 80/85/90% confidence
%                 ellipsoids per cluster and per reference class
%                 (DRAWGAUSSIANELLIPSOID), in the same colors as their
%                 own points. Skipped (silently) for any cluster/class
%                 with fewer than 4 points, same as PLOTPCARESULTS/
%                 PLOTGMMRESULTS's own 3D ellipsoids.
%     'Axes'      target axes (default: creates a new figure/axes). Must
%                 already be 3D-capable -- a fresh one created here
%                 always is.
%
%   Note: IDENTBESTDIST itself was computed in the full retained PCA
%   subspace (which can have more than 3 dimensions), while this plot
%   only has 3 axes to draw in -- a line's length here is a 3D projection
%   of that distance, not the distance itself.
%
%   Errors if the file has no reference-library identification result
%   (predates the feature, no reference library was loaded/overlaid for
%   that run, or fewer than 3 principal components were retained).
%
%   Example:
%     plotIdentification3DResults('run1.mat');
%     plotIdentification3DResults('run1.mat', 'Ellipses', false);
%
%   See also LOADPCARESULTS, IDENTIFYCLUSTERSBYREFERENCE,
%   PLOTREFERENCEOVERLAY, PLOTPCARESULTS, PCACLUSTERAPP.
    p = inputParser;
    addParameter(p, 'Ellipses', true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
    addParameter(p, 'Axes', [], @(v) isempty(v) || (isscalar(v) && isgraphics(v, 'axes')));
    parse(p, varargin{:});
    opt = p.Results;

    S = loadPCAResults(source);
    if isempty(S.refScore) || isempty(S.identBestClass)
        error('plotIdentification3DResults:noIdentification', ...
            ['This file has no reference-library identification result (it predates the ' ...
             'feature, or no reference library was loaded/overlaid for that run).']);
    end
    if size(S.score, 2) < 3
        error('plotIdentification3DResults:tooFewPCs', ...
            'This file only retained %d principal component(s); need at least 3 for a 3D view.', ...
            size(S.score, 2));
    end

    k = max(S.clusterIdx);
    if isempty(opt.Axes)
        fig = figure('Position', [100 100 900 700]);
        ax = axes(fig);
    else
        ax = opt.Axes;
    end
    cla(ax);

    clusterColors = lines(k);
    view(ax, 3);
    grid(ax, 'on');
    hold(ax, 'on');
    for c = 1:k
        mask = S.clusterIdx == c;
        scatter3(ax, S.score(mask,1), S.score(mask,2), S.score(mask,3), 36, clusterColors(c,:), ...
            'filled', 'DisplayName', sprintf('Cluster %d', c));
    end
    plotReferenceOverlay(ax, S.refScore(:,1:3), S.refClassUsed);

    if opt.Ellipses
        for c = 1:k
            mask = S.clusterIdx == c;
            if nnz(mask) >= 4
                drawGaussianEllipsoid(ax, mean(S.score(mask,1:3), 1), cov(S.score(mask,1:3)), clusterColors(c,:));
            end
        end
        refClassNames = unique(S.refClassUsed, 'stable');
        refColors = hsv(numel(refClassNames));  % same order/palette PLOTREFERENCEOVERLAY uses, so colors match
        for rc = 1:numel(refClassNames)
            mask = strcmp(S.refClassUsed, refClassNames{rc});
            if nnz(mask) >= 4
                drawGaussianEllipsoid(ax, mean(S.refScore(mask,1:3), 1), cov(S.refScore(mask,1:3)), refColors(rc,:));
            end
        end
    end

    for c = 1:k
        clusterCentroid = mean(S.score(S.clusterIdx == c, 1:3), 1);
        refMask = strcmp(S.refClassUsed, S.identBestClass{c});
        refCentroid = mean(S.refScore(refMask, 1:3), 1);
        plot3(ax, [clusterCentroid(1) refCentroid(1)], [clusterCentroid(2) refCentroid(2)], ...
            [clusterCentroid(3) refCentroid(3)], '--', 'Color', clusterColors(c,:), 'LineWidth', 2, ...
            'DisplayName', sprintf('Cluster %d -> %s (d=%.3f)', c, S.identBestClass{c}, S.identBestDist(c)));
        plot3(ax, clusterCentroid(1), clusterCentroid(2), clusterCentroid(3), 'o', ...
            'MarkerSize', 10, 'MarkerFaceColor', clusterColors(c,:), 'MarkerEdgeColor', 'k', ...
            'HandleVisibility', 'off');
        plot3(ax, refCentroid(1), refCentroid(2), refCentroid(3), 'd', ...
            'MarkerSize', 10, 'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerEdgeColor', 'k', ...
            'HandleVisibility', 'off');
    end
    hold(ax, 'off');
    legend(ax, 'Location', 'bestoutside');
    xlabel(ax, sprintf('PC1 (%.1f%%)', S.explained(1)));
    ylabel(ax, sprintf('PC2 (%.1f%%)', S.explained(2)));
    zlabel(ax, sprintf('PC3 (%.1f%%)', S.explained(3)));
    ax.XLimMode = 'auto';
    ax.YLimMode = 'auto';
    ax.ZLimMode = 'auto';
    title(ax, 'Cluster centroids vs. nearest reference class (PC1-PC2-PC3)');
end
