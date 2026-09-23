function ari = adjustedRandIndex(labelsA, labelsB)
%ADJUSTEDRANDINDEX  Agreement between two categorical partitions of the
%   same items, corrected for chance (Hubert & Arabie, 1985).
%
%   ARI = ADJUSTEDRANDINDEX(LABELSA, LABELSB) takes two label vectors of
%   the same length (numeric, categorical, or cell array of char/string,
%   in any encoding -- only which items share a label matters, not the
%   label values themselves) and returns a scalar in roughly [-1, 1]:
%   1 for identical partitions, ~0 for agreement no better than random,
%   negative for worse than random. Used here to check how well
%   unsupervised k-means clusters recover the filename-derived groups,
%   without that comparison influencing the clustering itself.
    [~, ~, a] = unique(labelsA(:));
    [~, ~, b] = unique(labelsB(:));
    n = numel(a);
    contingency = accumarray([a, b], 1);

    nchoose2 = @(x) x .* (x - 1) / 2;
    sumIJ = sum(nchoose2(contingency(:)));
    sumA = sum(nchoose2(sum(contingency, 2)));
    sumB = sum(nchoose2(sum(contingency, 1)));
    expectedIndex = sumA * sumB / nchoose2(n);
    maxIndex = 0.5 * (sumA + sumB);

    if maxIndex == expectedIndex
        ari = 0;  % every item in its own singleton partition on both sides
    else
        ari = (sumIJ - expectedIndex) / (maxIndex - expectedIndex);
    end
end
