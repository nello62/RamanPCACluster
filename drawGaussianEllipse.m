function drawGaussianEllipse(ax, mu, C, color)
%DRAWGAUSSIANELLIPSE  Draw 80/85/90% confidence ellipses of a bivariate
%   normal distribution.
%
%   DRAWGAUSSIANELLIPSE(AX, MU, C, COLOR) draws, into axes AX, the
%   80% (dotted), 85% (dashed), and 90% (solid) confidence ellipses of a
%   bivariate normal with mean MU (2-element vector) and covariance C
%   (2x2 matrix), in COLOR -- a standard way to visualize how tight or
%   overlapping a group of 2D points is, beyond the raw scatter alone.
%   MU/C can come from a group's own empirical mean/covariance (e.g. a
%   filename-derived group, or a k-means cluster), or directly from a
%   fitted model's own parameters (e.g. one GMM component) -- either way
%   the ellipse shows that specific distribution's shape, not a
%   recomputation from whichever points happen to be nearby.
%
%   Drawn with HandleVisibility off, so these never add clutter entries
%   to a legend already populated by GSCATTER/PLOT calls on the same axes
%   -- callers that want the ellipse to appear in a legend anyway are
%   expected to add their own entry (e.g. a dummy line with the desired
%   DisplayName) rather than rely on this one.
%
%   See also DRAWGAUSSIANELLIPSOID, PLOTPCARESULTS, PCACLUSTERAPP.
    [V, D] = eig(C);
    theta = linspace(0, 2*pi, 100);
    circle = [cos(theta); sin(theta)];
    confLevels = [0.80 0.85 0.90];
    styles = {':', '--', '-'};
    for i = 1:numel(confLevels)
        r = sqrt(chi2inv(confLevels(i), 2));
        pts = mu(:) + V * sqrt(D) * r * circle;
        plot(ax, pts(1,:), pts(2,:), styles{i}, 'Color', color, ...
            'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end
