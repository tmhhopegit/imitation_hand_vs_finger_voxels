function V = to_volume(values, voxels, dims, fill)
%TO_VOLUME Put per-voxel values back into a 3-D volume.
%   V = to_volume(values, voxels, dims)          untested voxels are NaN
%   V = to_volume(values, voxels, dims, fill)    ... or fill
if nargin < 4, fill = NaN; end
V = fill * ones(dims);
V(voxels) = double(values);
end
