function R = vlsm_analysis(D, cfg)
%VLSM_ANALYSIS Voxel-based lesion-symptom maps for two tasks and their difference.
%   R = vlsm_analysis(D, cfg)
%   D from load_imitation_data (or legacy_cell_to_data); cfg from config_voxel_maps.
%
%   Steps (replacing ImiHandFinger_StructureFunctionAnalysis, BatesVLSM,
%   BayesVLSM, WXVLSM, LPD_LPDD / LPD_Rudrauf2008):
%   1. Keep patients with both scores (and all covariates, if used).
%   2. Normalise the scores (normalise_scores: covariates regressed out,
%      rank-based inverse normal), and form the difference task1 - task2.
%   3. Keep voxels lesioned in at least cfg.vlsm.min_lesion_fraction of the
%      patients (and at least cfg.vlsm.min_lesion_count).
%   4. For each method in cfg.vlsm.methods and each target (task 1, task 2,
%      difference): the map, FDR threshold (only tested voxels count) and,
%      if cfg.vlsm.n_perm > 0, a family-wise threshold from permuting whole
%      patients.
%
%   R.maps is a struct array with fields target, method, stat (V x 1 over
%   R.voxels), p (for tests with p-values), fdr_threshold, fdr_mask,
%   fwer_threshold, fwer_mask. R.voxels are linear indices into D.dims;
%   use save_vlsm_maps to write them as images.
o = cfg.vlsm;
names = [D.score_names(:)', {sprintf('%s minus %s', D.score_names{1}, D.score_names{2})}];
cov = covariate_matrix(D, cfg.scores.covariates);
ok = all(~isnan(D.scores(:, 1:2)), 2);
if ~isempty(cov), ok = ok & ~any(isnan(cov), 2); end
if cfg.verbose, fprintf('%d of %d patients have both scores%s.\n', nnz(ok), numel(ok), ...
        conditional(~isempty(cov), ' and all covariates', '')); end
lesions = D.lesions(ok, :);
raw = D.scores(ok, 1:2);
if ~isempty(cov), cov = cov(ok, :); end
Y = normalise_scores(raw, cov, cfg.scores.method);
switch o.difference
    case 'of_normalised', Y(:, 3) = Y(:, 1) - Y(:, 2);                  % as Achilles.m
    case 'normalised_difference', Y(:, 3) = normalise_scores(raw(:, 1) - raw(:, 2), cov, cfg.scores.method);   % as Achilles2017.m
    otherwise, error('imitation:difference', 'cfg.vlsm.difference must be ''of_normalised'' or ''normalised_difference''.');
end
N = nnz(ok);
counts = sum(lesions, 1);
keep = counts >= max(o.min_lesion_count, ceil(o.min_lesion_fraction * N)) & counts <= N - o.min_lesion_count;
lesions = lesions(:, keep);
if cfg.verbose, fprintf('%d voxels are lesioned in enough patients to test.\n', nnz(keep)); end

R.patients = find(ok);
R.voxels = D.voxels(keep);
R.dims = D.dims;
R.scores = Y;
R.target_names = names;
R.maps = struct('target', {}, 'method', {}, 'stat', {}, 'p', {}, 'fdr_threshold', {}, ...
    'fdr_mask', {}, 'fwer_threshold', {}, 'fwer_mask', {});
for m = 1:numel(o.methods)
    method = o.methods{m};
    switch method
        case 'ttest'
            f = @(Yp) stats_for(@(y) vlsm_ttest(lesions, y, o.ttest_vartype), Yp);
            [S, P] = f(Y);
        case 'ranksum'
            f = @(Yp) stats_for(@(y) vlsm_ranksum(lesions, y), Yp);
            [S, P] = f(Y);
        case 'bayes'
            S = [vlsm_bayes(lesions, Y(:, 1), o.bayes_prior), vlsm_bayes(lesions, Y(:, 2), o.bayes_prior), ...
                 vlsm_bayes(lesions, Y(:, 3), o.bayes_prior)];
            P = nan(size(S));
            f = [];
        case 'lpd'
            imp = raw < cfg.scores.impairment_cutoffs(:)';        % impaired = below cut-off
            f = @(I) [vlsm_lpd(lesions, I(:, 1)), vlsm_lpd(lesions, I(:, 2)), ...
                      vlsm_lpd(lesions, I(:, 1)) - vlsm_lpd(lesions, I(:, 2))];
            S = f(imp);
            P = nan(size(S));
        otherwise
            error('imitation:method', 'Unknown VLSM method "%s".', method);
    end
    fwer = nan(1, 3);
    if o.n_perm > 0 && ~isempty(f)
        if cfg.verbose, fprintf('%s: %d permutations...\n', method, o.n_perm); end
        if strcmp(method, 'lpd'), perm_in = imp; else, perm_in = Y; end
        fwer = permutation_threshold(@(Yp) first_output(f, Yp), perm_in, o.n_perm, o.alpha, cfg.seed, cfg.workers);
    end
    for t = 1:3
        k = numel(R.maps) + 1;
        R.maps(k).target = names{t};
        R.maps(k).method = method;
        R.maps(k).stat = S(:, t);
        R.maps(k).p = P(:, t);
        if all(isnan(P(:, t)))
            R.maps(k).fdr_threshold = NaN;
            R.maps(k).fdr_mask = false(size(P(:, t)));
        else
            [R.maps(k).fdr_threshold, R.maps(k).fdr_mask] = fdr_threshold(P(:, t), o.fdr_q, o.fdr_method);
        end
        R.maps(k).fwer_threshold = fwer(t);
        R.maps(k).fwer_mask = abs(S(:, t)) >= fwer(t);
    end
end
end

function [S, P] = stats_for(fun, Y)
S = nan(0, size(Y, 2)); P = S;
for t = 1:size(Y, 2)
    [s, p] = fun(Y(:, t));
    S(1:numel(s), t) = s;
    P(1:numel(p), t) = p;
end
end

function S = first_output(f, Y)
S = f(Y);
end

function cov = covariate_matrix(D, names)
cov = [];
if isempty(names), return; end
[found, idx] = ismember(names, D.covariate_names);
if ~all(found)
    error('imitation:covariate', 'Unknown covariate(s): %s', strjoin(names(~found), ', '));
end
cov = D.covariates(:, idx);
end

function s = conditional(c, a, b)
if c, s = a; else, s = b; end
end
