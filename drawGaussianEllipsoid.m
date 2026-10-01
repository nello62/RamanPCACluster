function drawGaussianEllipsoid(ax, mu, C, color)
%DRAWGAUSSIANELLIPSOID  Draw 80/85/90% confidence ellipsoids of a
%   trivariate normal distribution.
%
%   DRAWGAUSSIANELLIPSOID(AX, MU, C, COLOR) draws, into 3D axes AX, three
%   nested semi-transparent confidence-ellipsoid shells (80/85/90%) of a
%   trivariate normal with mean MU (3-element vector) and covariance C
%   (3x3 matrix), in COLOR -- the 3D analogue of DRAWGAUSSIANELLIPSE's 2D
%   confidence ellipses, for the PC1-PC2-PC3 view.
%
%   Since three overlapping opaque surfaces would simply hide the
%   innermost one (unlike 2D's dotted/dashed/solid outlines, which stay
%   visible regardless of draw order), the levels are drawn largest-first
%   (90% down to 80%) with increasing opacity, so the tighter, more
%   confident region reads as the visually denser "core" seen through the
%   looser outer shell rather than being hidden by it.
%
%   HandleVisibility is off, so these never add legend entries.
%
%   See also DRAWGAUSSIANELLIPSE, PLOTPCARESULTS, PCACLUSTERAPP.
    [V, D] = eig(C);
    [xs, ys, zs] = sphere(30);
    unitSphere = [xs(:)'; ys(:)'; zs(:)'];

    confLevels = [0.90 0.85 0.80];  % largest (most transparent) first
    faceAlphas = [0.08 0.14 0.22];
    for i = 1:numel(confLevels)
        r = sqrt(chi2inv(confLevels(i), 3));
        pts = mu(:) + V * sqrt(D) * r * unitSphere;
        X = reshape(pts(1,:), size(xs));
        Y = reshape(pts(2,:), size(ys));
        Z = reshape(pts(3,:), size(zs));
        surf(ax, X, Y, Z, 'FaceColor', color, 'FaceAlpha', faceAlphas(i), ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end
