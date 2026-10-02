function [refWN, refIntensity, refClass, refNames] = loadReferenceSpectra(refDir)
%LOADREFERENCESPECTRA  Load a folder of reference Raman spectra (known
%   material standards, e.g. the SLoPP/SLoPP-E plastics library) for
%   visual comparison against a PCA analysis -- never for the analysis
%   itself.
%
%   [REFWN, REFINTENSITY, REFCLASS, REFNAMES] = LOADREFERENCESPECTRA(REFDIR)
%   reads every *.txt, *.dpt, *.csv, *.spc, and *.wdf file directly inside
%   REFDIR (non-recursive: a library folder organized both as a flat list
%   and as per-material subfolders of the same files, like SLoPP, would
%   otherwise double-count every spectrum), dispatched by extension to
%   the matching reader -- the same readers LOADRAMANSPECTRA and
%   RamanFitApp.m use:
%     .txt/.csv  two columns (wavenumber, intensity), no header
%                (delimiter auto-detected by READMATRIX)
%     .dpt       comma-separated wavenumber,intensity pairs, no header
%                (READDPT)
%     .spc       Galactic/GRAMS binary spectrum format (READSPC)
%     .wdf       Renishaw WiRE binary spectrum format (READWDF)
%   Unlike LOADRAMANSPECTRA's analysis spectra, reference spectra are NOT
%   required to share a common wavenumber grid (a library assembled from
%   many instruments/sessions usually won't, and mixing formats within
%   one library only makes that more likely), so each is kept on its own
%   native grid regardless of file format.
%
%   Returns, one entry per spectrum:
%     REFWN         - {n x 1} cell array, each cell the spectrum's own
%                     wavenumber vector
%     REFINTENSITY  - {n x 1} cell array, each cell the matching
%                     intensity vector
%     REFCLASS      - {n x 1} cellstr, the material class guessed from
%                     the file name: everything before the first
%                     " <number>." run (e.g. "Acrylic 1. Red Fiber
%                     (Polyacrylonitrile).txt" -> "Acrylic"), so that
%                     multiple reference spectra of the same material
%                     share one class label for grouping/coloring.
%     REFNAMES      - {n x 1} cellstr, the full file name (without
%                     extension) for per-spectrum labels/tooltips.
%
%   Non-spectrum files that happen to sit in REFDIR (e.g. SLoPP's own
%   "lib_names*.txt" index files, which are one name per line, not
%   wavenumber/intensity pairs) are skipped rather than raising an error:
%   any file whose name starts with "lib_" is skipped outright, and any
%   remaining file that fails to read as its extension's expected format,
%   or does not parse as an [N x 2] numeric matrix (N >= 2) for the
%   .txt/.csv case, is skipped with a warning -- so the function stays
%   robust to whatever else a library folder happens to contain.
%
%   See also PROJECTREFERENCESPECTRA, LOADRAMANSPECTRA, PCACLUSTERAPP.
    exts = {'*.txt', '*.dpt', '*.csv', '*.spc', '*.wdf'};
    files = [];
    for i = 1:numel(exts)
        files = [files; dir(fullfile(refDir, exts{i}))]; %#ok<AGROW>
    end
    filenames = {files.name}';
    filenames(startsWith(filenames, 'lib_', 'IgnoreCase', true)) = [];
    if isempty(filenames)
        error('loadReferenceSpectra:noFiles', ...
            'No .txt/.dpt/.csv/.spc/.wdf reference files found in %s.', refDir);
    end
    [~, ord] = sort(filenames);
    filenames = filenames(ord);
    n = numel(filenames);

    refWN = cell(n, 1);
    refIntensity = cell(n, 1);
    refClass = cell(n, 1);
    refNames = cell(n, 1);
    keep = true(n, 1);

    for i = 1:n
        f = fullfile(refDir, filenames{i});
        [~, ~, ext] = fileparts(f);
        try
            switch lower(ext)
                case '.dpt'
                    [x, y] = readdpt(f);
                case '.spc'
                    [x, y] = readspc(f);
                case '.wdf'
                    [x, y] = readwdf(f);
                otherwise  % .txt, .csv
                    data = readmatrix(f);
                    if ~isnumeric(data) || size(data, 2) ~= 2 || size(data, 1) < 2
                        error('loadReferenceSpectra:badShape', 'not a two-column spectrum');
                    end
                    x = data(:, 1);
                    y = data(:, 2);
            end
        catch
            warning('loadReferenceSpectra:skipped', ...
                '%s does not look like a two-column spectrum -- skipped.', filenames{i});
            keep(i) = false;
            continue
        end
        [wn, ordxy] = sort(x(:));
        refWN{i} = wn;
        refIntensity{i} = y(ordxy);
        [~, name, ~] = fileparts(filenames{i});
        refNames{i} = name;
        refClass{i} = regexprep(name, '^(.*?)\s+[0-9].*$', '$1');
    end

    refWN = refWN(keep);
    refIntensity = refIntensity(keep);
    refClass = refClass(keep);
    refNames = refNames(keep);
end
