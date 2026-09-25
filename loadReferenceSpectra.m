function [refWN, refIntensity, refClass, refNames] = loadReferenceSpectra(refDir)
%LOADREFERENCESPECTRA  Load a folder of reference Raman spectra (known
%   material standards, e.g. the SLoPP/SLoPP-E plastics library) for
%   visual comparison against a PCA analysis -- never for the analysis
%   itself.
%
%   [REFWN, REFINTENSITY, REFCLASS, REFNAMES] = LOADREFERENCESPECTRA(REFDIR)
%   reads every *.txt file directly inside REFDIR (non-recursive: a
%   library folder organized both as a flat list and as per-material
%   subfolders of the same files, like SLoPP, would otherwise double-count
%   every spectrum). Each file is two whitespace-separated columns
%   (wavenumber, intensity), no header -- unlike LOADRAMANSPECTRA's
%   .dpt files, reference spectra are NOT required to share a common
%   wavenumber grid (a library assembled from many instruments/sessions
%   usually won't), so each is kept on its own native grid.
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
%   Non-spectrum text files that happen to sit in REFDIR (e.g. SLoPP's own
%   own "lib_names*.txt" index files, which are one name per line, not
%   wavenumber/intensity pairs) are skipped rather than raising an error:
%   any file whose name starts with "lib_" is skipped outright, and any
%   remaining file that does not parse as an [N x 2] numeric matrix
%   (N >= 2) is skipped with a warning -- so the function stays robust to
%   whatever else a library folder happens to contain.
%
%   See also PROJECTREFERENCESPECTRA, LOADRAMANSPECTRA, PCACLUSTERAPP.
    files = dir(fullfile(refDir, '*.txt'));
    filenames = {files.name}';
    filenames(startsWith(filenames, 'lib_', 'IgnoreCase', true)) = [];
    if isempty(filenames)
        error('loadReferenceSpectra:noFiles', 'No .txt reference files found in %s.', refDir);
    end
    n = numel(filenames);

    refWN = cell(n, 1);
    refIntensity = cell(n, 1);
    refClass = cell(n, 1);
    refNames = cell(n, 1);
    keep = true(n, 1);

    for i = 1:n
        try
            data = readmatrix(fullfile(refDir, filenames{i}));
        catch
            data = [];
        end
        if ~isnumeric(data) || size(data, 2) ~= 2 || size(data, 1) < 2
            warning('loadReferenceSpectra:skipped', ...
                '%s does not look like a two-column spectrum -- skipped.', filenames{i});
            keep(i) = false;
            continue
        end
        [wn, ord] = sort(data(:, 1));
        refWN{i} = wn;
        refIntensity{i} = data(ord, 2);
        [~, name, ~] = fileparts(filenames{i});
        refNames{i} = name;
        refClass{i} = regexprep(name, '^(.*?)\s+[0-9].*$', '$1');
    end

    refWN = refWN(keep);
    refIntensity = refIntensity(keep);
    refClass = refClass(keep);
    refNames = refNames(keep);
end
