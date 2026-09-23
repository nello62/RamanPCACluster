function [Xproc, baselines] = preprocessSpectra(X, wavenumbers, options)
%PREPROCESSSPECTRA  Baseline-correct, smooth, and normalize a matrix of
%   Raman spectra before PCA/clustering.
%
%   XPROC = PREPROCESSSPECTRA(X, WAVENUMBERS) applies, to each row of X
%   (one spectrum per row):
%     1. Fluorescence baseline removal (one of backcor/airPLS/SNIP/APLS
%        -- the same four methods RamanFitApp offers), so PCA isn't
%        dominated by broad background differences between measurements
%        rather than real Raman peaks.
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
%     .Baseline        (default true)     set false to skip baseline removal
%     .BaselineMethod  (default 'airPLS') 'backcor' | 'airPLS' | 'SNIP' | 'APLS'
%     .BaselineLambda  (default 1e6)      airPLS smoothness parameter
%     .BaselineOrder   (default 2)        airPLS difference-penalty order
%     .BackcorOrder    (default 5)        backcor polynomial order
%     .BackcorThreshold (default 0.1)     backcor cost-function threshold
%     .BackcorFct      (default 'atq')    backcor cost function: sh|ah|stq|atq
%     .SnipIter        (default 40)       SNIP iterations/max clipping distance
%     .SnipUseLLS      (default true)     SNIP log-log-sqrt transform
%     .AplsGamma       (default 1e5)      APLS smoothness parameter
%     .AplsOrder       (default 2)        APLS difference-penalty order (1 or 2)
%     .AplsIter        (default 10)       APLS max reweighting iterations
%     .SmoothWindow    (default 9)        Savitzky-Golay window (odd)
%     .SmoothOrder     (default 3)        Savitzky-Golay polynomial order
%     .Smooth          (default true)     set false to skip smoothing
%     .Normalize       (default 'area')   'area' | 'snv' | 'max' | 'none'
%
%   See also BACKCOR, AIRPLS, SNIP, APLS, SGOLAYFILT, LOADRAMANSPECTRA.
    if nargin < 3
        options = struct();
    end
    opt = mergeOptions(options, struct( ...
        'Baseline', true, ...
        'BaselineMethod', 'airPLS', ...
        'BaselineLambda', 1e6, ...
        'BaselineOrder', 2, ...
        'BackcorOrder', 5, ...
        'BackcorThreshold', 0.1, ...
        'BackcorFct', 'atq', ...
        'SnipIter', 40, ...
        'SnipUseLLS', true, ...
        'AplsGamma', 1e5, ...
        'AplsOrder', 2, ...
        'AplsIter', 10, ...
        'SmoothWindow', 9, ...
        'SmoothOrder', 3, ...
        'Smooth', true, ...
        'Normalize', 'area'));

    [nSpectra, nPoints] = size(X);

    if opt.Baseline
        switch lower(opt.BaselineMethod)
            case 'airpls'
                % AIRPLS accepts the whole (nSpectra x nPoints) matrix
                % directly, one spectrum per row, and loops internally --
                % no need to call it per row like the other three methods.
                [Xproc, baselines] = airPLS(X, opt.BaselineLambda, opt.BaselineOrder);
            case 'backcor'
                baselines = zeros(nSpectra, nPoints);
                for i = 1:nSpectra
                    baselines(i, :) = backcor(wavenumbers, X(i, :), ...
                        opt.BackcorOrder, opt.BackcorThreshold, opt.BackcorFct)';
                end
                Xproc = X - baselines;
            case 'snip'
                baselines = zeros(nSpectra, nPoints);
                for i = 1:nSpectra
                    y = X(i, :);
                    if opt.SnipUseLLS && min(y) <= -1
                        % The LLS ("improved SNIP") transform computes
                        % log(log(sqrt(y+1)+1)+1), which goes complex for
                        % y <= -1; real detector data can dip slightly
                        % negative (dark-current offset), so a spectrum
                        % that needs it gets shifted to be safely
                        % non-negative first, with the same shift
                        % subtracted back out of the resulting baseline
                        % so it still lines up with the ORIGINAL spectrum.
                        shift = -min(y) + 1;
                        baselines(i, :) = snip(y + shift, opt.SnipIter, true) - shift;
                    else
                        baselines(i, :) = snip(y, opt.SnipIter, opt.SnipUseLLS);
                    end
                end
                Xproc = X - baselines;
            case 'apls'
                baselines = zeros(nSpectra, nPoints);
                for i = 1:nSpectra
                    baselines(i, :) = apls(X(i, :), opt.AplsGamma, opt.AplsOrder, opt.AplsIter);
                end
                Xproc = X - baselines;
            otherwise
                error('preprocessSpectra:badBaselineMethod', ...
                    'Unknown BaselineMethod: %s', opt.BaselineMethod);
        end
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
