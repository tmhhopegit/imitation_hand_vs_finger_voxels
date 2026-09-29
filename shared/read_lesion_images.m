function [lesions, voxels, dims, template] = read_lesion_images(files, opts)
%READ_LESION_IMAGES Binary lesion images as a patients x voxels matrix.
%   [lesions, voxels, dims, template] = read_lesion_images(files, opts)
%   opts.brain_mask  image whose non-zero voxels are kept ('' = all voxels)
%   opts.preprocess  options for preprocess_lesion, or [] for none
%   lesions: N x V logical; voxels: linear indices into a volume of size dims;
%   template: SPM header of the first image (for writing maps).
%   Every image must be in the same space (e.g. 2 mm MNI, 91 x 109 x 91).
%   Needs SPM.
if nargin < 2, opts = struct(); end
template = spm_vol(files{1});
dims = template.dim;
if isfield(opts, 'brain_mask') && ~isempty(opts.brain_mask)
    M = spm_read_vols(spm_vol(opts.brain_mask));
    if ~isequal(size(M), dims), error('imitation:space', 'The brain mask is not in the lesions'' space.'); end
    voxels = find(M > 0);
else
    voxels = (1:prod(dims))';
end
pre = [];
if isfield(opts, 'preprocess'), pre = opts.preprocess; end
lesions = false(numel(files), numel(voxels));
for i = 1:numel(files)
    V = spm_vol(files{i});
    if ~isequal(V.dim, dims), error('imitation:space', '%s is not in the same space as %s.', files{i}, files{1}); end
    X = spm_read_vols(V) > 0;
    if ~isempty(pre), X = preprocess_lesion(X, pre); end
    lesions(i, :) = X(voxels)';
end
end
