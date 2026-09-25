function [ldaScores, explainedLDA, cvAccuracy, confMat, classNames] = runLDA(X, groupLabels)
%RUNLDA  Supervised Linear Discriminant Analysis: canonical discriminant
%   scores for visualization, plus cross-validated classification
%   accuracy against known group labels.
%
%   [LDASCORES, EXPLAINEDLDA, CVACCURACY, CONFMAT, CLASSNAMES] =
%   RUNLDA(X, GROUPLABELS) takes X ([nSamples x nFeatures], typically
%   already-reduced PCA scores rather than raw spectra -- LDA's
%   within-class scatter matrix is singular whenever there are more
%   features than samples, which raw spectra almost always are) and one
%   group label per row of X (cell array of char/string, or categorical),
%   and returns:
%     LDASCORES    - [nSamples x min(nClasses-1, nFeatures)] canonical
%                    discriminant scores (LD1, LD2, ...), analogous to
%                    PCA's SCORE but supervised: these maximize
%                    separation between the given groups specifically,
%                    rather than total variance.
%     EXPLAINEDLDA - percentage of between-class variance captured by
%                    each discriminant axis (like PCA's EXPLAINED).
%     CVACCURACY   - k-fold cross-validated classification accuracy
%                    (fraction correctly classified on held-out folds
%                    the classifier never trained on) -- a measure of
%                    how separable the groups actually are, not just how
%                    they happen to look on a 2D projection of them.
%     CONFMAT      - confusion matrix from that same cross-validation.
%     CLASSNAMES   - class labels in CONFMAT's row/column order and
%                    LDASCORES' underlying group assignment order.
%
%   Requires at least two distinct groups in GROUPLABELS, and folds
%   capped at the smallest group's own sample count (a group with only 2
%   members can only support a 2-fold split, etc.) so no fold ends up
%   missing an entire class.
%
%   See also FITCDISCR, PCA, LOADRAMANSPECTRA.
    groupLabels = categorical(groupLabels(:));
    classNames = categories(groupLabels);
    nClasses = numel(classNames);
    if nClasses < 2
        error('runLDA:tooFewGroups', 'Need at least two distinct groups for LDA (found %d).', nClasses);
    end

    % --- Cross-validated classification accuracy ---
    mdl = fitcdiscr(X, groupLabels);
    kfold = max(2, min([10; countcats(groupLabels)]));
    cvmdl = crossval(mdl, 'KFold', kfold);
    cvAccuracy = 1 - kfoldLoss(cvmdl);
    predicted = kfoldPredict(cvmdl);
    confMat = confusionmat(groupLabels, predicted, 'Order', categorical(classNames));

    % --- Canonical discriminant scores (classical multi-class LDA) ---
    nFeatures = size(X, 2);
    overallMean = mean(X, 1);
    Sw = zeros(nFeatures);
    Sb = zeros(nFeatures);
    for c = 1:nClasses
        mask = groupLabels == classNames{c};
        Xc = X(mask, :);
        nc = size(Xc, 1);
        muC = mean(Xc, 1);
        Xc0 = Xc - muC;
        Sw = Sw + Xc0' * Xc0;
        d = (muC - overallMean)';
        Sb = Sb + nc * (d * d');
    end
    nComponents = min(nClasses - 1, nFeatures);
    [V, D] = eig(Sb, Sw);
    eigvals = diag(D);
    [eigvals, ord] = sort(eigvals, 'descend');
    V = V(:, ord(1:nComponents));
    eigvals = eigvals(1:nComponents);
    eigvals(eigvals < 0) = 0;  % numerical noise on near-zero eigenvalues can go slightly negative
    if sum(eigvals) > 0
        explainedLDA = 100 * eigvals / sum(eigvals);
    else
        explainedLDA = zeros(nComponents, 1);
    end
    ldaScores = (X - overallMean) * V;
end
