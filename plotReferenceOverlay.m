function plotReferenceOverlay(ax, refScore, refClassUsed)
%PLOTREFERENCEOVERLAY  Overlay already-projected reference spectra on a
%   PCA scatter plot, one marker shape per material class.
%
%   PLOTREFERENCEOVERLAY(AX, REFSCORE, REFCLASSUSED) draws REFSCORE(:,1:2)
%   (from PROJECTREFERENCESPECTRA -- reference spectra fitted INTO an
%   already-fitted PCA space, never part of fitting it) into axes AX as
%   small, distinctly-marked points, cycling through a fixed marker set
%   combined with a dedicated HSV color per class in REFCLASSUSED, so they
%   read as "known standards for comparison" rather than being mistaken
%   for actual data points/clusters. Call AX's HOLD ON before, and restore
%   its XLimMode/YLimMode to 'auto' after (a prior GSCATTER on the same
%   axes fixes them to just its own points, which will not automatically
%   grow to fit these) -- left to the caller since both are shared with
%   whatever else (e.g. confidence ellipses) the caller may also be
%   overlaying in the same HOLD block.
%
%   HandleVisibility stays on (the default) so each class gets its own
%   legend entry -- the entire point of the overlay -- unlike a purely
%   decorative addition to the plot.
%
%   See also PROJECTREFERENCESPECTRA, LOADREFERENCESPECTRA, GSCATTER.
    markers = {'o','s','d','^','v','>','<','p','h','x','+','*'};
    refClassNames = unique(refClassUsed, 'stable');
    refColors = hsv(numel(refClassNames));
    for c = 1:numel(refClassNames)
        m = markers{mod(c - 1, numel(markers)) + 1};
        mask = strcmp(refClassUsed, refClassNames{c});
        plot(ax, refScore(mask,1), refScore(mask,2), m, ...
            'Color', refColors(c,:), 'LineStyle', 'none', 'MarkerSize', 7, ...
            'LineWidth', 1.2, 'DisplayName', ['Ref: ' refClassNames{c}]);
    end
end
