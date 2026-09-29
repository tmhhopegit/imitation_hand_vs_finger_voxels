function C = compare_covariate_control(D, cfg)
%COMPARE_COVARIATE_CONTROL How much do the maps depend on which covariates are controlled?
%   C = compare_covariate_control(D, cfg)
%   Re-runs the voxel maps once per covariate set in cfg.covariate_sets
%   (rows of {name, {covariate names}}), e.g. no control, the covariates
%   used in the paper, and the full set. This is the Achilles analysis
%   (Achilles.m, Achilles2017.m, BatesVLSM_Raw_vs_Controlled, GetThreshMaps).
%
%   C.results{s}   vlsm_analysis output for set s
%   C.summary      table: set, target, method, voxels tested, FDR threshold,
%                  significant voxels, mean statistic of significant voxels
%   C.overlap      table comparing each set with the first: voxels
%                  significant only in the first, only in the other, in both,
%                  and the correlation between the two maps
%
%   All sets use the same patients (those with every covariate of every
%   set) and the same voxels, so the maps differ only by the control.
sets = cfg.covariate_sets;
allcov = unique([sets{:, 2}]);
if ~isempty(allcov)
    [found, idx] = ismember(allcov, D.covariate_names);
    if ~all(found), error('imitation:covariate', 'Unknown covariate(s): %s', strjoin(allcov(~found), ', ')); end
    complete = ~any(isnan(D.covariates(:, idx)), 2);
    D = subset_patients(D, complete);
end
C.results = cell(size(sets, 1), 1);
rows = {};
for s = 1:size(sets, 1)
    c = cfg;
    c.scores.covariates = sets{s, 2};
    if cfg.verbose, fprintf('--- covariate set "%s"\n', sets{s, 1}); end
    C.results{s} = vlsm_analysis(D, c);
    for m = C.results{s}.maps
        sig = m.fdr_mask;
        rows(end + 1, :) = {sets{s, 1}, m.target, m.method, numel(m.stat), m.fdr_threshold, nnz(sig), ...
            mean(m.stat(sig))}; %#ok<AGROW>
    end
end
C.summary = cell2table(rows, 'VariableNames', {'Set', 'Target', 'Method', 'VoxelsTested', ...
    'FDRThreshold', 'Significant', 'MeanStatSignificant'});

rows = {};
base = C.results{1};
for s = 2:size(sets, 1)
    other = C.results{s};
    if ~isequal(base.voxels, other.voxels)
        warning('imitation:voxels', 'Sets tested different voxels; overlaps use the common ones.');
    end
    [~, ia, ib] = intersect(base.voxels, other.voxels);
    for k = 1:numel(base.maps)
        a = base.maps(k); b = other.maps(k);
        ma = a.fdr_mask(ia); mb = b.fdr_mask(ib);
        rows(end + 1, :) = {sets{1, 1}, sets{s, 1}, a.target, a.method, nnz(ma & ~mb), nnz(~ma & mb), ...
            nnz(ma & mb), corr(a.stat(ia), b.stat(ib), 'rows', 'complete')}; %#ok<AGROW>
    end
end
C.overlap = cell2table(rows, 'VariableNames', {'SetA', 'SetB', 'Target', 'Method', 'OnlyA', 'OnlyB', ...
    'Both', 'MapCorrelation'});
end

function D = subset_patients(D, keep)
fields = {'ids', 'scores', 'covariates', 'lesions', 'loads', 'lesion_files'};
for f = fields
    if isfield(D, f{1}) && size(D.(f{1}), 1) == numel(keep)
        D.(f{1}) = D.(f{1})(keep, :);
    end
end
end
