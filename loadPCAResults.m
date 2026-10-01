function S = loadPCAResults(source)
%LOADPCARESULTS  Load (or pass through) a RamanPCACluster results struct,
%   normalized so every optional field is present.
%
%   S = LOADPCARESULTS(SOURCE) where SOURCE is a .mat file path (as saved
%   by PCAClusterApp.m's "Save results..." or by PCA_kmeans_analysis.m)
%   loads it with LOAD. SOURCE may also already be a loaded struct (e.g.
%   S0 = load('run1.mat'); S = loadPCAResults(S0)), returned as-is apart
%   from the same normalization -- every PLOT* function in this project
%   built on top of this one accepts either a file path or an in-memory
%   struct interchangeably.
%
%   A results file saved before the GMM/reference-library/identification
%   features existed won't have those fields at all; they are filled in
%   here as empty, so callers can check e.g. ISEMPTY(S.refScore) rather
%   than ISFIELD(S,'refScore') everywhere.
%
%   See also PLOTPCARESULTS, PCACLUSTERAPP, PCA_KMEANS_ANALYSIS.
    if ischar(source) || isstring(source)
        S = load(char(source));
    elseif isstruct(source)
        S = source;
    else
        error('loadPCAResults:badInput', 'SOURCE must be a .mat file path or a struct.');
    end

    required = {'score', 'explained', 'labels', 'clusterIdx'};
    missing = required(~isfield(S, required));
    if ~isempty(missing)
        error('loadPCAResults:missingFields', ...
            'Not a recognized RamanPCACluster results file/struct -- missing: %s.', ...
            strjoin(missing, ', '));
    end

    optionalDefaults = struct( ...
        'scoreReduced', [], ...
        'gmmClusterIdx', [], ...
        'refScore', [], 'refClassUsed', {{}}, 'refNamesUsed', {{}}, 'refSkipped', {{}}, ...
        'identBestClass', {{}}, 'identBestDist', [], ...
        'identSecondClass', {{}}, 'identSecondDist', []);
    fn = fieldnames(optionalDefaults);
    for i = 1:numel(fn)
        if ~isfield(S, fn{i})
            S.(fn{i}) = optionalDefaults.(fn{i});
        end
    end
end
