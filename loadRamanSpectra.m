function [X, wavenumbers, labels, filenames] = loadRamanSpectra(spectraDir)
%LOADRAMANSPECTRA  Load every .dpt Raman spectrum in a directory onto a
%   common wavenumber grid.
%
%   [X, WAVENUMBERS, LABELS, FILENAMES] = LOADRAMANSPECTRA(SPECTRADIR)
%   reads every *.dpt file in SPECTRADIR (READDPT format: comma-separated
%   wavenumber,intensity pairs, no header), returning:
%     X           - [nSpectra x nPoints] intensity matrix, one row per
%                   spectrum, sorted by increasing wavenumber
%     wavenumbers - [1 x nPoints] common wavenumber axis
%     labels      - {nSpectra x 1} group name guessed from each file's
%                   name (the part before the trailing "_<number><letter>.
%                   <number>.dpt" suffix, e.g. "CT_P_11b.0.dpt" -> "CT_P");
%                   purely for later validation/plotting, not used to
%                   guide the (unsupervised) analysis itself
%     filenames   - {nSpectra x 1} original file names, same order as X
%
%   All files must share the exact same wavenumber grid as the first one
%   read (sorted, compared to within 1e-6): this analysis stacks spectra
%   directly into a matrix rather than interpolating, so a file recorded
%   on a different grid raises an error instead of silently corrupting
%   the comparison.
    files = dir(fullfile(spectraDir, '*.dpt'));
    if isempty(files)
        error('loadRamanSpectra:noFiles', 'No .dpt files found in %s.', spectraDir);
    end
    filenames = {files.name}';
    nSpectra = numel(filenames);

    [x0, y0] = readdpt(fullfile(spectraDir, filenames{1}));
    [wavenumbers, ord0] = sort(x0(:)');
    nPoints = numel(wavenumbers);
    X = zeros(nSpectra, nPoints);
    y0 = y0(:)';
    X(1, :) = y0(ord0);

    for i = 2:nSpectra
        [xi, yi] = readdpt(fullfile(spectraDir, filenames{i}));
        [xi, ordi] = sort(xi(:)');
        if numel(xi) ~= nPoints || max(abs(xi - wavenumbers)) > 1e-6
            error('loadRamanSpectra:gridMismatch', ...
                '%s has a different wavenumber grid than %s.', filenames{i}, filenames{1});
        end
        yi = yi(:)';
        X(i, :) = yi(ordi);
    end

    labels = regexprep(filenames, '_[0-9]+[a-z]*\.[0-9]+\.dpt$', '');
end
