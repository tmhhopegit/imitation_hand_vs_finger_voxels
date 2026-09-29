function out = normalise_scores(scores, covariates, method)
%NORMALISE_SCORES Covariate-adjusted, normally distributed scores (from Norm_Resid_MinusThresh).
%   out = normalise_scores(scores, covariates, method)
%   scores: N x S; covariates: N x C (or [] for none).
%   method: 'rank_normal' (default): z-score, regress out the z-scored
%             covariates, then map the residuals' ranks onto a standard
%             normal (rank / (n + 1), then norminv). This is what the
%             original did.
%           'residual': the residuals without the rank step.
%           'none': the scores unchanged.
%   Patients with a missing score or covariate get NaN (the original's
%   zscore/regress returned NaN for everyone if any value was missing).
%   The original's third argument (imp_thresh) was never used.
if nargin < 2, covariates = []; end
if nargin < 3, method = 'rank_normal'; end
out = nan(size(scores));
if strcmp(method, 'none'), out = scores; return; end
for s = 1:size(scores, 2)
    ok = ~isnan(scores(:, s));
    if ~isempty(covariates), ok = ok & ~any(isnan(covariates), 2); end
    y = zscore(scores(ok, s));
    if ~isempty(covariates)
        X = [ones(nnz(ok), 1), zscore(covariates(ok, :))];
        y = y - X * (X \ y);
    end
    if strcmp(method, 'rank_normal')
        y = norminv(tiedrank(y) / (numel(y) + 1));
    elseif ~strcmp(method, 'residual')
        error('imitation:method', 'method must be ''rank_normal'', ''residual'' or ''none''.');
    end
    out(ok, s) = y;
end
end
