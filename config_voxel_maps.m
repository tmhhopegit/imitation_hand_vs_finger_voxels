function cfg = config_voxel_maps()
%CONFIG_VOXEL_MAPS Settings for the voxel-based lesion-symptom maps.
%   cfg = config_voxel_maps();  then change fields as needed.
%   Everything the original scripts hard-coded (column numbers, image size,
%   cut-offs, FDR level, permutations) is set here.

% ---------------------------------------------------------------- data (load_imitation_data)
cfg.data.spreadsheet        = '';                    % one row per patient, with a header row
cfg.data.id_column          = 'ID';
cfg.data.score_columns      = {'Hand', 'Finger'};    % the two imitation scores (task 1, task 2)
cfg.data.covariate_columns  = {'Age', 'TimePost', 'LesionVolume'};   % 'Name=A+B' = sum of two columns, called Name
cfg.data.lesion_folder      = '';                    % binary lesion images, searched recursively
cfg.data.lesion_pattern     = '*.nii';
cfg.data.lesion_id_regex    = '^r?(?<id>[^_]+)';     % patient ID = file name up to the first '_' (after an optional 'r')
cfg.data.brain_mask         = '';                    % keep only these voxels ('' = all)
cfg.data.preprocess         = [];                    % e.g. struct('gap_fill_slices', 6, 'one_hemisphere', true) as GetLesions
cfg.data.region_folder      = '';                    % region images for lesion loads ('' = none)
cfg.data.region_pattern     = '*.nii';
cfg.data.load_columns       = {};                    % or region loads from the spreadsheet: names, or the first such column

% ---------------------------------------------------------------- scores
cfg.scores.covariates         = {};              % covariates regressed out of both scores (names from covariate_columns)
cfg.scores.method             = 'rank_normal';   % 'rank_normal' (Norm_Resid_MinusThresh), 'residual' or 'none'
cfg.scores.impairment_cutoffs = [18 17];         % impaired = score below this (hand, finger); for 'lpd' maps

% ---------------------------------------------------------------- maps
cfg.vlsm.methods             = {'ttest', 'bayes'};   % any of 'ttest', 'ranksum', 'bayes', 'lpd'
cfg.vlsm.min_lesion_fraction = 0.1;              % test voxels lesioned in at least 10% of patients (MinSample)
cfg.vlsm.min_lesion_count    = 2;                % ... and in at least this many, and not in all but this many
cfg.vlsm.difference          = 'of_normalised';  % difference map: task1 - task2 of the normalised scores (Achilles.m)
                                                 % or 'normalised_difference': normalise the raw difference (Achilles2017.m)
cfg.vlsm.ttest_vartype       = 'unequal';        % Welch, as BatesVLSM/Achilles; 'equal' as WXVLSM's copy
cfg.vlsm.bayes_prior         = 'jzs';            % 'jzs' or 'unit_information' (bayes_anova1's prior, corrected)
cfg.vlsm.fdr_q               = 0.01;             % FDR level (the Achilles scripts used 1%)
cfg.vlsm.fdr_method          = 'by';             % 'by' = fdr_bh 'dep' (any dependence), 'bh' = 'pdep'
cfg.vlsm.n_perm              = 0;                % permutations for family-wise thresholds (e.g. 1000; 0 = none)
cfg.vlsm.alpha               = 0.05;             % family-wise error rate

% ---------------------------------------------------------------- covariate comparison (Achilles)
% rows: {name, covariates}; the first set is the reference for the overlaps
cfg.covariate_sets = {
    'no control',  {}
    'as in paper', {'Age', 'TimePost'}
    'full',        {'Age', 'TimePost', 'LesionVolume'}};

% ---------------------------------------------------------------- general
cfg.seed        = 1;
cfg.workers     = 0;          % parfor workers for permutations (0 = serial)
cfg.verbose     = true;
cfg.results_dir = fullfile(pwd, 'results');
end
