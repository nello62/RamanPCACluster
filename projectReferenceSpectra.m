function [refScore, refClassUsed, refNamesUsed, skipped] = projectReferenceSpectra( ...
    refWN, refIntensity, refClass, refNames, wn, opt, meanXproc, coeff)
%PROJECTREFERENCESPECTRA  Map reference spectra into an already-fitted PCA
%   space, for visual comparison only -- they never take part in fitting
%   PCA/k-means/GMM/LDA, which is the whole point of a reference library
%   (known material standards to compare the actual clustering against,
%   not additional data to cluster).
%
%   [REFSCORE, REFCLASSUSED, REFNAMESUSED, SKIPPED] =
%   PROJECTREFERENCESPECTRA(REFWN, REFINTENSITY, REFCLASS, REFNAMES, WN,
%   OPT, MEANXPROC, COEFF) takes the per-spectrum cell arrays from
%   LOADREFERENCESPECTRA (each reference on its own native wavenumber
%   grid) and, for every one that covers the analysis range:
%     1. Skips it if its own native grid does not fully cover
%        [min(WN) max(WN)] -- extrapolating a Raman baseline or peak
%        beyond what was actually measured for that reference is not
%        reliable, so it is left out rather than silently guessed at.
%     2. Interpolates it (linear) onto WN, the exact grid the main
%        analysis was run on.
%     3. Runs it through PREPROCESSSPECTRA with OPT, the SAME options
%        struct used for the main dataset, so baseline removal/smoothing/
%        normalization are done identically -- otherwise "the same PCA
%        space" would not mean the same thing for the two data sets.
%     4. Projects it via (Xproc - MEANXPROC) * COEFF, reusing the
%        loadings COEFF and the centering mean MEANXPROC from the
%        already-fitted PCA (i.e. MEAN(Xproc,1) of the analyzed dataset,
%        not of the references) -- consistent with how PCA's own SCORE
%        output is defined, but never refitting PCA on the references.
%
%   Returns:
%     REFSCORE      - [nUsed x nComponents] projected scores, nComponents
%                     = size(COEFF,2) (as many columns as COEFF has)
%     REFCLASSUSED, REFNAMESUSED - the class/name of each kept reference,
%                     same order/length as REFSCORE's rows
%     SKIPPED       - {k x 1} cellstr, names of references left out for
%                     not covering the analysis range (informational, so
%                     the caller can report it rather than the reference
%                     just silently vanishing)
%
%   See also LOADREFERENCESPECTRA, PREPROCESSSPECTRA, PCA, PCACLUSTERAPP.
    n = numel(refWN);
    wnMin = min(wn);
    wnMax = max(wn);

    covers = false(n, 1);
    for i = 1:n
        covers(i) = min(refWN{i}) <= wnMin && max(refWN{i}) >= wnMax;
    end
    skipped = refNames(~covers);

    idxUsed = find(covers);
    nUsed = numel(idxUsed);
    Xref = zeros(nUsed, numel(wn));
    for j = 1:nUsed
        i = idxUsed(j);
        Xref(j, :) = interp1(refWN{i}, refIntensity{i}, wn, 'linear');
    end

    if nUsed == 0
        refScore = zeros(0, size(coeff, 2));
        refClassUsed = {};
        refNamesUsed = {};
        return
    end

    XprocRef = preprocessSpectra(Xref, wn, opt);
    refScore = (XprocRef - meanXproc) * coeff;
    refClassUsed = refClass(idxUsed);
    refNamesUsed = refNames(idxUsed);
end
