function [clusterIdx, post, bic, gmModel, bicScan, kScanUsed] = runGMM(X, k, kScan)
%RUNGMM  Gaussian Mixture Model soft clustering via Expectation-Maximization.
%
%   [CLUSTERIDX, POST, BIC, GMMODEL] = RUNGMM(X, K) fits a K-component
%   Gaussian Mixture Model to X ([nSamples x nFeatures], typically the
%   same PCA-reduced scores used for k-means) via FITGMDIST (full,
%   unshared covariance per component; regularized against a singular
%   covariance estimate for components that end up with few or nearly
%   collinear points). Returns:
%     CLUSTERIDX - [nSamples x 1] hard cluster assignment (argmax of the
%                  posterior), for a like-for-like comparison against
%                  k-means' own CLUSTERIDX.
%     POST       - [nSamples x K] posterior probability of each component
%                  given each point -- the "soft" part of soft
%                  clustering, unlike k-means' all-or-nothing assignment.
%     BIC        - Bayesian Information Criterion of the fitted model at
%                  K components (lower is better), from GMMODEL.BIC.
%     GMMODEL    - the fitted gmdistribution object (component means,
%                  covariances, mixing proportions), used by the caller
%                  to draw each component's covariance ellipse.
%
%   [..., BICSCAN, KSCANUSED] = RUNGMM(X, K, KSCAN) additionally refits
%   the model at every component count in KSCAN (default 1:8, capped to
%   at most NSAMPLES-1) and returns each fit's BIC in BICSCAN, aligned
%   with KSCANUSED -- a model-selection diagnostic analogous to k-means'
%   elbow/silhouette scan, except this one comes from the GMM's own
%   likelihood rather than a separate criterion. A component count that
%   fails to fit (e.g. a collapsing/near-singular covariance) is recorded
%   as NaN in BICSCAN rather than aborting the whole scan.
%
%   See also FITGMDIST, PCA, KMEANS.
    if nargin < 3
        kScan = 1:8;
    end
    nSamples = size(X, 1);
    kScan = kScan(kScan <= nSamples - 1 & kScan >= 1);

    regVal = 1e-6;
    opts = statset('MaxIter', 500);

    gmModel = fitgmdist(X, k, 'CovarianceType', 'full', ...
        'RegularizationValue', regVal, 'Replicates', 5, 'Options', opts);
    post = posterior(gmModel, X);
    [~, clusterIdx] = max(post, [], 2);
    bic = gmModel.BIC;

    if nargout > 4
        bicScan = nan(size(kScan));
        for ik = 1:numel(kScan)
            try
                mdlIk = fitgmdist(X, kScan(ik), 'CovarianceType', 'full', ...
                    'RegularizationValue', regVal, 'Replicates', 3, 'Options', opts);
                bicScan(ik) = mdlIk.BIC;
            catch
                bicScan(ik) = NaN;
            end
        end
        kScanUsed = kScan;
    end
end
