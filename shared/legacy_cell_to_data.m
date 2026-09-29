function D = legacy_cell_to_data(Data, cols, brain_vox, dims)
%LEGACY_CELL_TO_DATA Convert the old Data cell array (ConnectLesionsToData) to the new format.
%   D = legacy_cell_to_data(Data, cols, brain_vox, dims)
%   Data: cell array with a header row; one patient per row after it.
%   cols: struct of column numbers, as the old scripts used them:
%     cols.id         = 1
%     cols.scores     = [14 15]        (hand, finger)
%     cols.covariates = [16 17 18 19]  (any demographic columns)
%     cols.lesion     = 13             (3-D lesion images)
%     cols.loads      = 21:size(Data,2) (region loads; optional)
%   brain_vox: linear indices of the voxels to keep (e.g. from GetMasks);
%   dims: volume size, default [91 109 91].
%   D.template is empty, so give save_volume a template (D.template =
%   spm_vol('any_lesion.nii')) before writing maps.
if nargin < 4, dims = [91 109 91]; end
rows = 2:size(Data, 1);
D.ids = cellstr(string(Data(rows, cols.id)));
D.score_names = cellstr(string(Data(1, cols.scores)));
D.scores = cell_numbers(Data(rows, cols.scores));
D.covariate_names = cellstr(string(Data(1, cols.covariates)));
D.covariates = cell_numbers(Data(rows, cols.covariates));
D.voxels = brain_vox(:);
D.dims = dims;
D.template = [];
D.lesions = false(numel(rows), numel(brain_vox));
for i = 1:numel(rows)
    X = Data{rows(i), cols.lesion};
    D.lesions(i, :) = X(brain_vox)' > 0;
end
if isfield(cols, 'loads') && ~isempty(cols.loads)
    D.region_names = cellstr(string(Data(1, cols.loads)));
    D.loads = cell_numbers(Data(rows, cols.loads));
end
end

function X = cell_numbers(C)
% empty cells (missing scores) become NaN; cell2mat would fail on them
X = nan(size(C));
for k = 1:numel(C)
    if ~isempty(C{k}) && isnumeric(C{k}), X(k) = double(C{k}); end
end
end
