function save_volume(values, D, file)
%SAVE_VOLUME Write per-voxel values as a NIfTI/Analyze image (needs SPM).
%   save_volume(values, D, file) uses D.voxels, D.dims and D.template
%   (from load_imitation_data). Untested voxels are written as NaN.
%   Replaces SaveVolume, which was not in the folder.
V = D.template;
V.fname = file;
V.dt = [spm_type('float32'), spm_platform('bigend')];
V.pinfo = [1; 0; 0];
if isfield(V, 'private'), V = rmfield(V, 'private'); end
spm_write_vol(V, to_volume(double(values), D.voxels, D.dims));
end
