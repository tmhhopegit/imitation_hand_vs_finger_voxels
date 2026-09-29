function R = interaction_analysis(D, improvement, factors, cfg)
%INTERACTION_ANALYSIS Voxels where the lesion effect depends on a factor (MakeInteractionMap_Thresh).
%   R = interaction_analysis(D, improvement, factors, cfg)
%   improvement: N x 1 outcome (e.g. change in score); factors: N x F, the
%   first interacts with the lesion. Patients with missing values, or any
%   factor <= 0 (as the original), are left out.
%   R.p, R.coef: per tested voxel (R.voxels); R.fwer_threshold: the p below
%   which a voxel is significant after permutation (cfg.vlsm.n_perm);
%   R.fwer_mask.
%   The original tested every voxel lesioned in at least one patient; this
%   uses cfg.vlsm.min_lesion_fraction / min_lesion_count like the other maps.
ok = ~isnan(improvement(:)) & min(factors, [], 2) > 0 & all(~isnan(factors), 2);
L = D.lesions(ok, :);
y = improvement(ok);
F = factors(ok, :);
N = nnz(ok);
counts = sum(L, 1);
keep = counts >= max(cfg.vlsm.min_lesion_count, ceil(cfg.vlsm.min_lesion_fraction * N)) & counts <= N - cfg.vlsm.min_lesion_count;
L = L(:, keep);
[R.p, R.coef] = vlsm_interaction(L, y, F);
R.voxels = D.voxels(keep);
R.dims = D.dims;
R.fwer_threshold = NaN;
R.fwer_mask = false(size(R.p));
if cfg.vlsm.n_perm > 0
    % permute the outcome only (factors stay with their patients), as the original
    [thr, ~] = permutation_threshold(@(yp) -log10(vlsm_interaction(L, yp, F)), y, ...
        cfg.vlsm.n_perm, cfg.vlsm.alpha, cfg.seed, cfg.workers);
    R.fwer_threshold = 10 ^ (-thr);
    R.fwer_mask = R.p <= R.fwer_threshold;
end
end
