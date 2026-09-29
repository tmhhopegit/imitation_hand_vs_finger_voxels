%EXAMPLE_VOXEL_MAPS Voxel-based lesion-symptom maps, section by section.
% Edit the paths and column names below, then run each section (Ctrl+Enter).
run(fullfile(fileparts(mfilename('fullpath')), 'setup_voxel_maps.m'));

cfg = config_voxel_maps();
cfg.data.spreadsheet   = 'C:\path\to\imitation_data.xlsx';
cfg.data.lesion_folder = 'C:\path\to\lesions';
cfg.data.brain_mask    = 'C:\path\to\brain_mask.nii';
cfg.data.covariate_columns = {'Age', 'TimePost=TimeYears+TimeMonths', 'LesionVolume'};   % as Achilles.m: d(:,4)+d(:,5)

%% Load the data (needs SPM)
D = load_imitation_data(cfg);
% Or from an old Data cell array (ConnectLesionsToData layout):
%   cols = struct('id', 1, 'scores', [14 15], 'covariates', 16:19, 'lesion', 13);
%   [~, bv] = GetMasks;  D = legacy_cell_to_data(Data, cols, bv);
%   D.template = spm_vol('C:\path\to\any_lesion.nii');

%% Maps for hand, finger and hand - finger, with FDR thresholds
cfg.scores.covariates = {'Age', 'TimePost', 'LesionVolume'};
cfg.vlsm.methods = {'ttest', 'ranksum', 'bayes', 'lpd'};
cfg.vlsm.n_perm = 1000;                  % family-wise thresholds too (slow-ish)
R = vlsm_analysis(D, cfg);
struct2table(rmfield(R.maps, {'stat', 'p', 'fdr_mask', 'fwer_mask'}))
save_vlsm_maps(R, D, fullfile(cfg.results_dir, 'maps'));

%% How much does the choice of covariates change the maps? (Achilles)
cfg.vlsm.methods = {'ttest'};
cfg.vlsm.n_perm = 0;
C = compare_covariate_control(D, cfg);
disp(C.summary)
disp(C.overlap)
