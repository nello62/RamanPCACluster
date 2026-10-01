function ax = plotScreeResults(source, varargin)
%PLOTSCREERESULTS  Re-plot the PCA scree plot (explained + cumulative
%   explained variance) from a saved RamanPCACluster results file,
%   without rerunning the analysis.
%
%   AX = PLOTSCREERESULTS(SOURCE) loads SOURCE (a .mat file path, or an
%   already-loaded struct -- see LOADPCARESULTS) and plots, into a new
%   figure, each principal component's own explained variance (bars, left
%   axis) alongside the cumulative explained variance (line, right axis)
%   -- the same plot as PCAClusterApp's own "PCA" tab, reconstructed from
%   the saved file alone.
%
%   AX = PLOTSCREERESULTS(..., 'Name', Value, ...) options:
%     'NumComponents'  how many leading PCs to show (default: min(10,
%                      total available) -- same default as the app).
%     'Axes'           target axes (default: creates a new figure/axes).
%
%   Examples:
%     plotScreeResults('run1.mat');
%     plotScreeResults('run1.mat', 'NumComponents', 15);
%
%   See also LOADPCARESULTS, PLOTPCARESULTS, PCACLUSTERAPP.
    p = inputParser;
    addParameter(p, 'NumComponents', [], @(v) isempty(v) || (isscalar(v) && isnumeric(v) && v >= 1));
    addParameter(p, 'Axes', [], @(v) isempty(v) || (isscalar(v) && isgraphics(v, 'axes')));
    parse(p, varargin{:});
    opt = p.Results;

    S = loadPCAResults(source);
    if isempty(opt.NumComponents)
        nShow = min(10, numel(S.explained));
    else
        nShow = min(round(opt.NumComponents), numel(S.explained));
    end

    if isempty(opt.Axes)
        fig = figure;
        ax = axes(fig);
    else
        ax = opt.Axes;
    end
    cla(ax);

    yyaxis(ax, 'left');
    bar(ax, S.explained(1:nShow));
    ylabel(ax, 'Explained variance (%)');
    yyaxis(ax, 'right');
    plot(ax, cumsum(S.explained(1:nShow)), '-o');
    ylabel(ax, 'Cumulative explained variance (%)');
    xlabel(ax, 'Principal component');
    title(ax, 'Scree plot');
end
