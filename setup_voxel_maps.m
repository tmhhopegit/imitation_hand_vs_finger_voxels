%SETUP_VOXEL_MAPS Put the voxel-map code on the MATLAB path.
%   Run once per session:  run('<this folder>/setup_voxel_maps.m')
%   The project is self-contained: its own copies of the shared helpers go
%   first on the path, so they are used even if the region project is on
%   the path too.
root = fileparts(mfilename('fullpath'));
for sub = {'tests', 'shared', 'analysis'}
    addpath(fullfile(root, sub{1}));
end
addpath(root);
clear root sub
