function [bestClass, bestDist, secondClass, secondDist] = identifyClustersByReference( ...
    clusterCentroids, refScoreSub, refClassUsed)
%IDENTIFYCLUSTERSBYREFERENCE  Indicative material identification for each
%   cluster, by nearest reference-class centroid in PCA space.
%
%   [BESTCLASS, BESTDIST, SECONDCLASS, SECONDDIST] =
%   IDENTIFYCLUSTERSBYREFERENCE(CLUSTERCENTROIDS, REFSCORESUB,
%   REFCLASSUSED) takes:
%     CLUSTERCENTROIDS - [nClusters x nDims] each cluster's own centroid
%                        (e.g. mean PCA score of its members), in
%                        whatever PCA subspace the caller chooses
%     REFSCORESUB      - [nRef x nDims] reference spectra already
%                        projected into that SAME subspace (see
%                        PROJECTREFERENCESPECTRA -- same number of
%                        columns/dimensions as CLUSTERCENTROIDS)
%     REFCLASSUSED     - {nRef x 1} cellstr, each reference's material
%                        class (see LOADREFERENCESPECTRA)
%
%   and returns, for every cluster:
%     BESTCLASS/BESTDIST     - the material class whose reference spectra
%                              (averaged into that class's own centroid)
%                              sit closest, and the Euclidean distance
%     SECONDCLASS/SECONDDIST - the runner-up class/distance, so the
%                              caller can gauge how confident or
%                              ambiguous the match is (e.g. flag it when
%                              the top two are nearly tied)
%
%   This is a purely descriptive heuristic (nearest-centroid matching, no
%   fitted decision boundary or held-out validation) meant as an
%   indicative starting point, not a validated classifier: "this cluster
%   sits closest, in PCA space, to these known reference spectra" -- a
%   small or incomplete reference library, or a cluster whose true
%   material simply is not represented in it, can still return a
%   "closest" class, which is then only the least-wrong guess among
%   what's in the library. Always read the reported distance (and the
%   PCA scatter itself) alongside the label, not the label alone.
%
%   See also PROJECTREFERENCESPECTRA, LOADREFERENCESPECTRA, PCACLUSTERAPP.
    refClassNames = unique(refClassUsed, 'stable');
    nClasses = numel(refClassNames);
    classCentroids = zeros(nClasses, size(refScoreSub, 2));
    for c = 1:nClasses
        classCentroids(c, :) = mean(refScoreSub(strcmp(refClassUsed, refClassNames{c}), :), 1);
    end

    nClusters = size(clusterCentroids, 1);
    bestClass = cell(nClusters, 1);
    bestDist = zeros(nClusters, 1);
    secondClass = cell(nClusters, 1);
    secondDist = nan(nClusters, 1);
    for i = 1:nClusters
        d = sqrt(sum((classCentroids - clusterCentroids(i, :)).^2, 2));
        [dSorted, ord] = sort(d);
        bestClass{i} = refClassNames{ord(1)};
        bestDist(i) = dSorted(1);
        if nClasses >= 2
            secondClass{i} = refClassNames{ord(2)};
            secondDist(i) = dSorted(2);
        else
            secondClass{i} = '';
        end
    end
end
