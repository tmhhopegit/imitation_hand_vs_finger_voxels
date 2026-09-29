function D = load_imitation_data(cfg)
%LOAD_IMITATION_DATA Scores, covariates, lesions and region loads for the analyses.
%   D = load_imitation_data(cfg) uses cfg.data:
%     spreadsheet        .xlsx/.csv, one row per patient, a header row
%     id_column          column with patient IDs
%     score_columns      the two task scores, e.g. {'Hand', 'Finger'}
%     covariate_columns  e.g. {'Age', 'TimePost=TimeYears+TimeMonths', 'LesionVolume'};
%                        'A+B' is the sum of columns A and B, and 'Name=...'
%                        names the result (otherwise the text itself is the name)
%     lesion_folder      folder of binary lesion images (searched recursively)
%     lesion_pattern     e.g. '*.nii'
%     lesion_id_regex    how to find the patient ID in a lesion file name;
%                        needs a group named id, e.g. '^r?(?<id>[^_]+)'
%     brain_mask         image defining the voxels to keep ('' = every voxel)
%     preprocess         options for preprocess_lesion (gap filling, folding
%                        onto one hemisphere); [] = none
%     region_folder      folder of region images for lesion loads ('' = none)
%     region_pattern     e.g. '*.nii'
%     load_columns       region loads already in the spreadsheet instead:
%                        a cell array of column names, or the name of the
%                        first load column (all columns from there on)
%
%   D has fields ids, scores (N x 2), score_names, covariates, covariate_names,
%   lesions (N x V logical), voxels (linear indices), dims, template (SPM
%   header for writing maps), lesion_files, loads (N x R), region_names,
%   region_masks (R x V logical).
%
%   Patients are matched to lesion files by EXACT ID (case-insensitive).
%   ConnectLesionsToData matched by substring, so ID 'p1' also matched
%   'p10', and an unmatched file crashed it. Unmatched files and patients
%   are listed in D.unmatched_files / D.missing_lesions.
%   Needs SPM (spm_vol, spm_read_vols) for images.
d = cfg.data;
T = readtable(d.spreadsheet, 'VariableNamingRule', 'preserve');
ids = cellstr(string(T.(d.id_column)));
D.ids = ids;
D.score_names = d.score_columns;
D.scores = column_values(T, d.score_columns);
[D.covariate_names, specs] = split_names(d.covariate_columns);
D.covariates = column_values(T, specs);

if isfield(d, 'lesion_folder') && ~isempty(d.lesion_folder)
    files = dir(fullfile(d.lesion_folder, '**', d.lesion_pattern));
    files = fullfile({files.folder}, {files.name})';
    file_ids = cell(size(files));
    for f = 1:numel(files)
        [~, name] = fileparts(files{f});
        tok = regexp(name, d.lesion_id_regex, 'names', 'once');
        if isempty(tok), file_ids{f} = ''; else, file_ids{f} = lower(tok.id); end
    end
    match = zeros(numel(ids), 1);
    for i = 1:numel(ids)
        k = find(strcmp(file_ids, lower(ids{i})));
        if numel(k) > 1
            error('imitation:ids', 'Patient %s matches %d lesion files.', ids{i}, numel(k));
        end
        if ~isempty(k), match(i) = k; end
    end
    D.missing_lesions = ids(match == 0);
    D.unmatched_files = files(~ismember(1:numel(files), match));
    if ~isempty(D.missing_lesions)
        warning('imitation:ids', '%d patients have no lesion image; they are dropped.', numel(D.missing_lesions));
    end
    keep = match > 0;
    D = select_rows(D, keep);
    D.lesion_files = files(match(keep));
    [D.lesions, D.voxels, D.dims, D.template] = read_lesion_images(D.lesion_files, d);
end

if isfield(d, 'region_folder') && ~isempty(d.region_folder)
    files = dir(fullfile(d.region_folder, '**', d.region_pattern));
    files = sort(fullfile({files.folder}, {files.name}))';
    [~, D.region_names] = cellfun(@fileparts, files, 'UniformOutput', false);
    D.region_masks = false(numel(files), numel(D.voxels));
    for r = 1:numel(files)
        R = spm_read_vols(spm_vol(files{r}));
        D.region_masks(r, :) = R(D.voxels)' > 0;
    end
    D.loads = region_loads(D.lesions, D.region_masks);
elseif isfield(d, 'load_columns') && ~isempty(d.load_columns)
    names = d.load_columns;
    if ischar(names)                     % first load column: take it and all after it
        vn = T.Properties.VariableNames;
        names = vn(find(strcmp(vn, names)):end);
    end
    D.region_names = names;
    D.loads = column_values(T, names);
    if isfield(D, 'rows'), D.loads = D.loads(D.rows, :); end
end
end

function [names, specs] = split_names(entries)
names = entries; specs = entries;
for c = 1:numel(entries)
    k = strfind(entries{c}, '=');
    if ~isempty(k)
        names{c} = strtrim(entries{c}(1:k(1) - 1));
        specs{c} = strtrim(entries{c}(k(1) + 1:end));
    end
end
end

function X = column_values(T, names)
X = nan(height(T), numel(names));
for c = 1:numel(names)
    parts = strtrim(strsplit(names{c}, '+'));
    for p = 1:numel(parts)
        if ~ismember(parts{p}, T.Properties.VariableNames)
            error('imitation:column', 'No column "%s" in the spreadsheet.', parts{p});
        end
        v = double(T.(parts{p}));
        if p == 1, X(:, c) = v; else, X(:, c) = X(:, c) + v; end
    end
end
end

function D = select_rows(D, keep)
D.rows = find(keep);
D.ids = D.ids(keep);
D.scores = D.scores(keep, :);
D.covariates = D.covariates(keep, :);
end
