function ax = plotPCAResults(source, varargin)
%PLOTPCARESULTS  Re-plot PCA scores from a saved RamanPCACluster results
%   file -- experimental data, confidence ellipses/ellipsoids, and any
%   saved reference-library projection -- without rerunning the analysis.
%
%   AX = PLOTPCARESULTS(SOURCE) loads SOURCE (a .mat file path, or an
%   already-loaded struct -- see LOADPCARESULTS) and plots PC1 vs PC2
%   into a new figure, colored by k-means cluster, with 80/85/90%
%   confidence ellipses per cluster and, if the file has one, the
%   reference-library overlay (PLOTREFERENCEOVERLAY) -- the same
%   information as PCAClusterApp's own "Clusters" tab, reproduced from
%   the saved file alone.
%
%   AX = PLOTPCARESULTS(..., 'Name', Value, ...) options:
%     'Dims'      [1 2] (default) or a 3-element vector, e.g. [1 2 3] --
%                 which principal components to plot. A 3-element DIMS
%                 switches to a 3D SCATTER3 plot (with confidence
%                 ellipsoids instead of ellipses) automatically.
%     'ColorBy'   'kmeans' (default) | 'group' | 'gmm' -- which partition
%                 colors the experimental points. Falls back to 'group'
%                 with a warning if the requested one isn't in the file
%                 (e.g. 'gmm' when the file predates GMM, or GMM failed
%                 to fit).
%     'Ellipses'  true (default) | false -- 80/85/90% confidence
%                 ellipses (2D) / ellipsoids (3D) per group/cluster.
%     'Reference' true (default) | false -- overlay the saved
%                 reference-library projection, if the file has one
%                 (silently skipped regardless of this value if it
%                 doesn't).
%     'Axes'      target axes handle (default: creates a new figure/axes).
%                 Must already be a 3D-capable axes if DIMS has 3
%                 elements -- a fresh one created here always is.
%
%   Examples:
%     plotPCAResults('run1.mat');                     % PC1 vs PC2, k-means
%     plotPCAResults('run1.mat', 'ColorBy', 'group');  % by filename group
%     plotPCAResults('run1.mat', 'ColorBy', 'gmm', 'Ellipses', false);
%     plotPCAResults('run1.mat', 'Dims', [1 2 3]);     % 3D PC1-PC2-PC3
%
%   See also LOADPCARESULTS, DRAWGAUSSIANELLIPSE, DRAWGAUSSIANELLIPSOID,
%   PLOTREFERENCEOVERLAY, PCACLUSTERAPP, PCA_KMEANS_ANALYSIS.
    p = inputParser;
    addParameter(p, 'Dims', [1 2], @(v) isnumeric(v) && (numel(v) == 2 || numel(v) == 3));
    addParameter(p, 'ColorBy', 'kmeans', @(v) any(strcmpi(v, {'kmeans', 'group', 'gmm'})));
    addParameter(p, 'Ellipses', true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
    addParameter(p, 'Reference', true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
    addParameter(p, 'Axes', [], @(v) isempty(v) || (isscalar(v) && isgraphics(v, 'axes')));
    parse(p, varargin{:});
    opt = p.Results;
    dims = round(opt.Dims(:)');
    is3D = numel(dims) == 3;

    S = loadPCAResults(source);
    if any(dims > size(S.score, 2))
        error('plotPCAResults:badDims', ...
            'This file only has %d principal components; requested Dims go up to %d.', ...
            size(S.score, 2), max(dims));
    end

    colorBy = lower(opt.ColorBy);
    if strcmp(colorBy, 'gmm') && isempty(S.gmmClusterIdx)
        warning('plotPCAResults:noGMM', 'No GMM result in this file -- falling back to ColorBy=''group''.');
        colorBy = 'group';
    end
    switch colorBy
        case 'kmeans'
            groupVar = S.clusterIdx;
            groupLabel = 'k-means cluster';
        case 'gmm'
            groupVar = S.gmmClusterIdx;
            groupLabel = 'GMM component';
        otherwise
            groupVar = S.labels;
            groupLabel = 'filename-derived group';
    end

    if isnumeric(groupVar)
        groupNames = num2cell(unique(groupVar));
    else
        groupNames = unique(groupVar);  % cellstr -> sorted alphabetically, same convention PCACLUSTERAPP uses
    end
    nGroups = numel(groupNames);
    colors = lines(nGroups);

    if isempty(opt.Axes)
        fig = figure;
        ax = axes(fig);
    else
        ax = opt.Axes;
    end
    cla(ax);
    hold(ax, 'on');

    axLabel = @(d) sprintf('PC%d (%.1f%%)', d, S.explained(d));

    if ~is3D
        gscatter(ax, S.score(:, dims(1)), S.score(:, dims(2)), groupVar, colors);
        if opt.Ellipses
            for g = 1:nGroups
                mask = matchesGroup(groupVar, groupNames{g});
                if nnz(mask) >= 3
                    mu = mean(S.score(mask, dims), 1);
                    C = cov(S.score(mask, dims(1)), S.score(mask, dims(2)));
                    drawGaussianEllipse(ax, mu, C, colors(g, :));
                end
            end
        end
        if opt.Reference && ~isempty(S.refScore)
            plotReferenceOverlay(ax, S.refScore(:, dims), S.refClassUsed);
        end
        xlabel(ax, axLabel(dims(1)));
        ylabel(ax, axLabel(dims(2)));
        ax.XLimMode = 'auto';
        ax.YLimMode = 'auto';
    else
        view(ax, 3);
        grid(ax, 'on');
        for g = 1:nGroups
            mask = matchesGroup(groupVar, groupNames{g});
            if isnumeric(groupVar)
                dispName = sprintf('%d', groupNames{g});
            else
                dispName = groupNames{g};
            end
            scatter3(ax, S.score(mask, dims(1)), S.score(mask, dims(2)), S.score(mask, dims(3)), ...
                36, colors(g, :), 'filled', 'DisplayName', dispName);
            if opt.Ellipses && nnz(mask) >= 4
                mu = mean(S.score(mask, dims), 1);
                C = cov(S.score(mask, dims));
                drawGaussianEllipsoid(ax, mu, C, colors(g, :));
            end
        end
        if opt.Reference && ~isempty(S.refScore)
            plotReferenceOverlay(ax, S.refScore(:, dims), S.refClassUsed);
        end
        legend(ax);
        xlabel(ax, axLabel(dims(1)));
        ylabel(ax, axLabel(dims(2)));
        zlabel(ax, axLabel(dims(3)));
        ax.XLimMode = 'auto';
        ax.YLimMode = 'auto';
        ax.ZLimMode = 'auto';
    end
    hold(ax, 'off');
    title(ax, sprintf('PC%s scores, colored by %s', strjoin(string(dims), '-'), groupLabel));
end

% -------------------------------------------------------------------------
function mask = matchesGroup(groupVar, name)
    if isnumeric(groupVar)
        mask = groupVar == name;
    else
        mask = strcmp(groupVar, name);
    end
end
