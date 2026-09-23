function [Xproc, baselines] = preprocessSpectra(X, wavenumbers, options)
%PREPROCESSSPECTRA  Baseline-correct, smooth, and normalize a matrix of
%   Raman spectra before PCA/clustering.
%
%   XPROC = PREPROCESSSPECTRA(X, WAVENUMBERS) applies, to each row of X
%   (one spectrum per row):
%     1. Fluorescence baseline removal (AIRPLS), so PCA isn't dominated
%        by broad background differences between measurements rather
%        than real Raman peaks.
%     2. Savitzky-Golay smoothing (light, order 3 / window 9 by default),
%        to reduce high-frequency detector noise that would otherwise
%        show up as spurious high-order PCs.
%     3. Row-wise normalization (unit area by default), so PCA isn't
%        dominated by overall signal-intensity differences between
%        measurements (laser power, focus, exposure) rather than
%        differences in Raman spectral shape.
%
%   [XPROC, BASELINES] = PREPROCESSSPECTRA(...) also returns the
%   estimated baseline for each spectrum (same size as X), for inspection.
%
%   OPTIONS is an optional struct overriding any of these fields:
%     .Baseline        (default true)  set false to skip baseline removal
%     .BaselineLambda  (default 1e6)   airPLS smoothness parameter
%     .BaselineOrder   (default 2)     airPLS difference order
%     .SmoothWindow    (default 9)     Savitzky-Golay window (odd)
%     .SmoothOrder     (default 3)     Savitzky-Golay polynomial order
%     .Smooth          (default true)  set false to skip smoothing
%     .Normalize       (default 'area') 'area' | 'snv' | 'max' | 'none'
%
%   See also AIRPLS, SGOLAYFILT, LOADRAMANSPECTRA.
    if nargin < 3
        options = struct();
    end
    opt = mergeOptions(options, struct( ...
        'Baseline', true, ...
        'BaselineLambda', 1e6, ...
        'BaselineOrder', 2, ...
        'SmoothWindow', 9, ...
        'SmoothOrder', 3, ...
        'Smooth', true, ...
        'Normalize', 'area'));

    [nSpectra, nPoints] = size(X);

    if opt.Baseline
        % AIRPLS accepts the whole (nSpectra x nPoints) matrix directly,
        % one spectrum per row, and loops internally -- no need to call
        % it per row.
        [Xproc, baselines] = airPLS(X, opt.BaselineLambda, opt.BaselineOrder);
    else
        Xproc = X;
        baselines = zeros(size(X));
    end

    if opt.Smooth
        for i = 1:nSpectra
            Xproc(i, :) = sgolayfilt(Xproc(i, :), opt.SmoothOrder, opt.SmoothWindow);
        end
    end

    switch lower(opt.Normalize)
        case 'area'
            % Negative dips left by baseline/noise would otherwise offset
            % the area; clipped at zero only for this normalizing
            % integral, not in the returned spectrum itself.
            areas = trapz(wavenumbers, max(Xproc, 0), 2);
            areas(areas <= 0) = 1;
            Xproc = Xproc ./ areas;
        case 'max'
            peaks = max(Xproc, [], 2);
            peaks(peaks <= 0) = 1;
            Xproc = Xproc ./ peaks;
        case 'snv'
            mu = mean(Xproc, 2);
            sigma = std(Xproc, 0, 2);
            sigma(sigma == 0) = 1;
            Xproc = (Xproc - mu) ./ sigma;
        case 'none'
            % leave as-is
        otherwise
            error('preprocessSpectra:badNormalize', 'Unknown Normalize option: %s', opt.Normalize);
    end
end

% -------------------------------------------------------------------------
function merged = mergeOptions(user, defaults)
    merged = defaults;
    f = fieldnames(user);
    for k = 1:numel(f)
        merged.(f{k}) = user.(f{k});
    end
end
