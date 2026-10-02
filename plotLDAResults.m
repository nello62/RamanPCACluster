function [axScatter, axConfusion] = plotLDAResults(source, varargin)
%PLOTLDARESULTS  Re-plot the LDA canonical discriminant scatter and
%   cross-validated confusion matrix from a saved RamanPCACluster results
%   file, without rerunning the analysis.
%
%   [AXSCATTER, AXCONFUSION] = PLOTLDARESULTS(SOURCE) loads SOURCE (a .mat
%   file path, or an already-loaded struct -- see LOADPCARESULTS) and
%   draws, into a new figure with two panels:
%     - LD1 vs LD2 (or LD1 vs. sample index if only one discriminant axis
%       exists, i.e. exactly two filename-derived groups), colored by
%       filename-derived group;
%     - the cross-validated confusion matrix (rows: actual group; columns:
%       predicted group), with per-cell counts and the overall
%       cross-validated accuracy in the title
%   -- the same information as PCAClusterApp's own "LDA" tab, reconstructed
%   from the saved file alone.
%
%   [...] = PLOTLDARESULTS(..., 'Name', Value, ...) options:
%     'Axes'  [axScatter, axConfusion], two target axes (default: creates
%             a new figure with a 1x2 tiled layout).
%
%   Errors if the file has no LDA result (it predates the feature, or
%   fewer than two filename-derived groups were present at analysis time --
%   ONRUNANALYSIS leaves LDASCORES empty in that case too).
%
%   Example:
%     plotLDAResults('run1.mat');
%
%   See also LOADPCARESULTS, PLOTPCARESULTS, PLOTGMMRESULTS, RUNLDA,
%   PCACLUSTERAPP.
    p = inputParser;
    addParameter(p, 'Axes', [], @(v) isempty(v) || (numel(v) == 2 && all(isgraphics(v, 'axes'))));
    parse(p, varargin{:});
    opt = p.Results;

    S = loadPCAResults(source);
    if isempty(S.ldaScores)
        error('plotLDAResults:noLDA', ...
            ['This file has no LDA result (it predates the feature, or fewer than two ' ...
             'filename-derived groups were present at analysis time).']);
    end

    if isempty(opt.Axes)
        fig = figure('Position', [100 100 1000 450]);
        tl = tiledlayout(fig, 1, 2);
        axScatter = nexttile(tl);
        axConfusion = nexttile(tl);
    else
        axScatter = opt.Axes(1);
        axConfusion = opt.Axes(2);
    end
    cla(axScatter);
    cla(axConfusion);

    groupColors = lines(numel(S.ldaClassNames));
    hold(axScatter, 'on');
    if size(S.ldaScores, 2) >= 2
        gscatter(axScatter, S.ldaScores(:,1), S.ldaScores(:,2), S.labels, groupColors);
        xlabel(axScatter, sprintf('LD1 (%.1f%%)', S.explainedLDA(1)));
        ylabel(axScatter, sprintf('LD2 (%.1f%%)', S.explainedLDA(2)));
    else
        % Exactly two groups: Fisher's LDA only has one discriminant axis,
        % so LD1 is plotted against sample index instead of a (nonexistent) LD2.
        gscatter(axScatter, (1:size(S.ldaScores,1))', S.ldaScores(:,1), S.labels, groupColors);
        xlabel(axScatter, 'Sample index');
        ylabel(axScatter, sprintf('LD1 (%.1f%%)', S.explainedLDA(1)));
    end
    hold(axScatter, 'off');
    axScatter.XLimMode = 'auto';
    axScatter.YLimMode = 'auto';
    title(axScatter, 'LD1 vs LD2 (colored by filename-derived group)');

    nC = numel(S.ldaClassNames);
    imagesc(axConfusion, S.ldaConfMat);
    colormap(axConfusion, 'parula');
    colorbar(axConfusion);
    axis(axConfusion, 'square');
    set(axConfusion, 'XTick', 1:nC, 'YTick', 1:nC, 'XTickLabel', S.ldaClassNames, ...
        'YTickLabel', S.ldaClassNames, 'XTickLabelRotation', 45, 'YDir', 'reverse');
    xlabel(axConfusion, 'Predicted (cross-validated)');
    ylabel(axConfusion, 'Actual (filename-derived group)');
    hold(axConfusion, 'on');
    for r = 1:nC
        for c = 1:nC
            text(axConfusion, c, r, num2str(S.ldaConfMat(r,c)), 'HorizontalAlignment', 'center', 'Color', [1 1 1]);
        end
    end
    hold(axConfusion, 'off');
    title(axConfusion, sprintf('Confusion matrix (cross-validated accuracy = %.1f%%)', 100 * S.cvAccuracy));
end
