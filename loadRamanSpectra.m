function [X, wavenumbers, labels, filenames] = loadRamanSpectra(spectraDir)
%LOADRAMANSPECTRA  Load every Raman spectrum in a directory onto a common
%   wavenumber grid.
%
%   [X, WAVENUMBERS, LABELS, FILENAMES] = LOADRAMANSPECTRA(SPECTRADIR)
%   reads every *.dpt, *.csv, *.spc, and *.wdf file in SPECTRADIR,
%   returning:
%     X           - [nSpectra x nPoints] intensity matrix, one row per
%                   spectrum, sorted by increasing wavenumber
%     wavenumbers - [1 x nPoints] common wavenumber axis
%     labels      - {nSpectra x 1} group name guessed from each file's
%                   name (everything up to the underscore that
%                   introduces the sample identifier, itself always
%                   starting with a digit -- e.g. "CT_P_11b.0.dpt" ->
%                   "CT_P", "Paradiso_31-nonpellet.0.dpt" -> "Paradiso";
%                   the sample identifier after that point can be
%                   anything, since real filenames in the wild are not
%                   always consistently formatted); purely for later
%                   validation/plotting, not used to guide the
%                   (unsupervised) analysis itself
%     filenames   - {nSpectra x 1} original file names, same order as X
%
%   Supported file formats, dispatched by extension (same readers as
%   RamanFitApp.m, so a folder exported for that app can be reused here
%   unchanged):
%     .dpt  comma-separated wavenumber,intensity pairs, no header
%           (READDPT; Hamamatsu/OPUS-style instrument export)
%     .csv  exactly two numeric columns (wavenumber, intensity), no
%           header row (plain READMATRIX; a header row or extra columns
%           are rejected rather than silently misread -- this also
%           catches a results .csv this app itself saved into the same
%           folder, which has a header and more columns)
%     .spc  Galactic/GRAMS binary spectrum format (READSPC, wrapping
%           GSSpcRead)
%     .wdf  Renishaw WiRE binary spectrum format (READWDF, wrapping
%           WdfReader)
%   Every file is expected to hold exactly one spectrum -- unlike
%   RamanFitApp, a "wide" multi-spectrum file (several intensity columns
%   sharing one wavenumber axis) is not supported here, since this
%   function assumes one file = one row of X = one filename-derived
%   label. Mixing formats within the same folder is fine as long as the
%   underlying wavenumber grid matches.
%
%   All files must share the exact same wavenumber grid as the first one
%   read (sorted, compared to within 1e-6): this analysis stacks spectra
%   directly into a matrix rather than interpolating, so a file recorded
%   on a different grid raises an error instead of silently corrupting
%   the comparison.
    exts = {'*.dpt', '*.csv', '*.spc', '*.wdf'};
    files = [];
    for i = 1:numel(exts)
        files = [files; dir(fullfile(spectraDir, exts{i}))]; %#ok<AGROW>
    end
    if isempty(files)
        error('loadRamanSpectra:noFiles', 'No .dpt/.csv/.spc/.wdf files found in %s.', spectraDir);
    end
    [~, ord] = sort({files.name});
    files = files(ord);
    filenames = {files.name}';
    nSpectra = numel(filenames);

    [x0, y0] = readSpectrumFile(fullfile(spectraDir, filenames{1}));
    [wavenumbers, ord0] = sort(x0(:)');
    nPoints = numel(wavenumbers);
    X = zeros(nSpectra, nPoints);
    y0 = y0(:)';
    X(1, :) = y0(ord0);

    for i = 2:nSpectra
        [xi, yi] = readSpectrumFile(fullfile(spectraDir, filenames{i}));
        [xi, ordi] = sort(xi(:)');
        if numel(xi) ~= nPoints || max(abs(xi - wavenumbers)) > 1e-6
            error('loadRamanSpectra:gridMismatch', ...
                '%s has a different wavenumber grid than %s.', filenames{i}, filenames{1});
        end
        yi = yi(:)';
        X(i, :) = yi(ordi);
    end

    % Extension stripped first (uniformly across .dpt/.csv/.spc/.wdf)
    % before the group-name pattern is applied, so a double extension
    % like the real ".0.dpt" naming doesn't need its own special case per
    % file format.
    [~, baseNames] = cellfun(@fileparts, filenames, 'UniformOutput', false);
    labels = regexprep(baseNames, '^(.*)_[0-9].*$', '$1');
end

% -------------------------------------------------------------------------
function [x, y] = readSpectrumFile(f)
    [~, ~, ext] = fileparts(f);
    switch lower(ext)
        case '.dpt'
            [x, y] = readdpt(f);
        case '.spc'
            [x, y] = readspc(f);
        case '.wdf'
            [x, y] = readwdf(f);
        case '.csv'
            data = readmatrix(f);
            % Deliberately strict (exactly two numeric columns, no
            % header row) rather than just ">= 2 columns": this is also
            % what catches a stray non-spectrum .csv sitting in the same
            % folder (e.g. this app's own *_cluster_assignments.csv,
            % *_reference_projections.csv output, which has 6+ columns
            % and a header row) with a clear, specific error instead of a
            % confusing downstream "different wavenumber grid" message or
            % silently misreading garbage columns as wavenumber/intensity.
            if size(data, 2) ~= 2 || any(isnan(data(:)))
                error('loadRamanSpectra:badFile', ...
                    ['Expected exactly two numeric columns (wavenumber, intensity), no header ' ...
                     'row, in %s -- got %d column(s) (or a non-numeric/header row). If this is ' ...
                     'not meant to be a spectrum (e.g. a previously saved results CSV), move it ' ...
                     'out of the spectra folder.'], f, size(data, 2));
            end
            x = data(:, 1)';
            y = data(:, 2)';
        otherwise
            error('loadRamanSpectra:badExtension', 'Unsupported file extension %s (%s).', ext, f);
    end
end
