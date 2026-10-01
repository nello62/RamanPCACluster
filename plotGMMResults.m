function [axScatter, axBIC] = plotGMMResults(source, varargin)
%PLOTGMMRESULTS  Re-plot the GMM soft-clustering scatter and BIC
%   model-selection scan from a saved RamanPCACluster results file,
%   without rerunning the analysis.
%
%   [AXSCATTER, AXBIC] = PLOTGMMRESULTS(SOURCE) loads SOURCE (a .mat file
%   path, or an already-loaded struct -- see LOADPCARESULTS) and draws,
%   into a new figure with two panels:
%     - PC1 vs PC2 (or 'Dims', see below), colored by GMM component
%       (hard, argmax-posterior assignment), with each component's own
%       confidence ellipse built from the FITTED MODEL's own mean/
%       covariance (GMMMODEL.mu/.Sigma) -- the model's actual shape, not
%       a post-hoc empirical summary of whichever points it ended up
%       being assigned;
%     - BIC vs. number of components (the saved scan across
%       k = 1,...,8), with the lowest-BIC point starred and the k
%       actually used marked by a dashed line
%   -- the same information as PCAClusterApp's own "GMM" tab,
%   reconstructed from the saved file alone.
%
%   [...] = PLOTGMMRESULTS(..., 'Name', Value, ...) options:
%     'Dims'      [1 2] (default) or a 3-element vector, e.g. [1 2 3] --
%                 which principal components to plot in the scatter
%                 panel. A 3-element DIMS switches it to a 3D SCATTER3
%                 view with confidence ellipsoids (see
%                 DRAWGAUSSIANELLIPSOID), same convention as
%                 PLOTPCARESULTS.
%     'Ellipses'  true (default) | false -- component confidence
%                 ellipses (2D) / ellipsoids (3D).
%     'Axes'      [axScatter, axBIC], two target axes (default: creates a
%                 new figure with a 1x2 tiled layout).
%
%   Errors if the file has no GMM result (e.g. it predates the feature,
%   or the fit failed for that run -- ONRUNANALYSIS leaves GMMCLUSTERIDX
%   empty in that case too).
%
%   Example:
%     plotGMMResults('run1.mat');
%     plotGMMResults('run1.mat', 'Dims', [1 2 3]);
%
%   See also LOADPCARESULTS, PLOTPCARESULTS, DRAWGAUSSIANELLIPSE,
%   DRAWGAUSSIANELLIPSOID, PCACLUSTERAPP.
    p = inputParser;
    addParameter(p, 'Dims', [1 2], @(v) isnumeric(v) && (numel(v) == 2 || numel(v) == 3));
    addParameter(p, 'Ellipses', true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
    addParameter(p, 'Axes', [], @(v) isempty(v) || (numel(v) == 2 && all(isgraphics(v, 'axes'))));
    parse(p, varargin{:});
    opt = p.Results;
    dims = round(opt.Dims(:)');
    is3D = numel(dims) == 3;

    S = loadPCAResults(source);
    if isempty(S.gmmClusterIdx)
        error('plotGMMResults:noGMM', ...
            'This file has no GMM result (it predates the feature, or the fit failed for that run).');
    end
    if any(dims > size(S.score, 2))
        error('plotGMMResults:badDims', ...
            'This file only has %d principal components; requested Dims go up to %d.', ...
            size(S.score, 2), max(dims));
    end

    nComp = size(S.gmmPosterior, 2);
    gmmColors = lines(nComp);
    meanMaxPost = mean(max(S.gmmPosterior, [], 2));

    if isempty(opt.Axes)
        fig = figure('Position', [100 100 1000 450]);
        tl = tiledlayout(fig, 1, 2);
        axScatter = nexttile(tl);
        axBIC = nexttile(tl);
    else
        axScatter = opt.Axes(1);
        axBIC = opt.Axes(2);
    end
    cla(axScatter);
    cla(axBIC);

    axLabel = @(d) sprintf('PC%d (%.1f%%)', d, S.explained(d));
    hold(axScatter, 'on');
    if ~is3D
        gscatter(axScatter, S.score(:, dims(1)), S.score(:, dims(2)), S.gmmClusterIdx, gmmColors);
        if opt.Ellipses
            for c = 1:nComp
                drawGaussianEllipse(axScatter, S.gmmModel.mu(c, dims), extractSigma(S.gmmModel, c, dims), gmmColors(c, :));
            end
        end
        xlabel(axScatter, axLabel(dims(1)));
        ylabel(axScatter, axLabel(dims(2)));
        axScatter.XLimMode = 'auto';
        axScatter.YLimMode = 'auto';
    else
        view(axScatter, 3);
        grid(axScatter, 'on');
        for c = 1:nComp
            mask = S.gmmClusterIdx == c;
            scatter3(axScatter, S.score(mask, dims(1)), S.score(mask, dims(2)), S.score(mask, dims(3)), ...
                36, gmmColors(c, :), 'filled', 'DisplayName', sprintf('%d', c));
            if opt.Ellipses
                drawGaussianEllipsoid(axScatter, S.gmmModel.mu(c, dims), extractSigma(S.gmmModel, c, dims), gmmColors(c, :));
            end
        end
        legend(axScatter);
        xlabel(axScatter, axLabel(dims(1)));
        ylabel(axScatter, axLabel(dims(2)));
        zlabel(axScatter, axLabel(dims(3)));
        axScatter.XLimMode = 'auto';
        axScatter.YLimMode = 'auto';
        axScatter.ZLimMode = 'auto';
    end
    hold(axScatter, 'off');
    title(axScatter, sprintf('GMM soft clustering (k=%d, mean max-posterior = %.2f)', nComp, meanMaxPost));

    plot(axBIC, S.gmmKScanUsed, S.gmmBICScan, '-o');
    hold(axBIC, 'on');
    [minBIC, ib] = min(S.gmmBICScan);
    if ~isnan(minBIC)
        plot(axBIC, S.gmmKScanUsed(ib), minBIC, 'r*', 'MarkerSize', 10, 'HandleVisibility', 'off');
    end
    yl = ylim(axBIC);
    plot(axBIC, [nComp nComp], yl, 'k--', 'HandleVisibility', 'off');
    hold(axBIC, 'off');
    xlabel(axBIC, 'Number of components');
    ylabel(axBIC, 'BIC (lower is better)');
    title(axBIC, sprintf('Model selection: BIC vs. components (used k=%d, dashed line)', nComp));
end

% -------------------------------------------------------------------------
function Sigma = extractSigma(gmmModel, c, dims)
    if strcmpi(gmmModel.CovarianceType, 'diagonal')
        Sigma = diag(gmmModel.Sigma(1, dims, c));
    else
        Sigma = gmmModel.Sigma(dims, dims, c);
    end
end
