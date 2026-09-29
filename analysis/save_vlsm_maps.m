function files = save_vlsm_maps(R, D, folder, what)
%SAVE_VLSM_MAPS Write the maps from vlsm_analysis as images (needs SPM).
%   files = save_vlsm_maps(R, D, folder)
%   files = save_vlsm_maps(R, D, folder, what)   what: any of 'stat', 'p',
%            'fdr', 'fwer' (default all). Thresholded maps hold the
%            statistic where significant and 0 elsewhere.
%   File names: <method>_<target>_<what>.nii
if nargin < 4, what = {'stat', 'p', 'fdr', 'fwer'}; end
if ~exist(folder, 'dir'), mkdir(folder); end
files = {};
T.voxels = R.voxels; T.dims = R.dims; T.template = D.template;
for m = R.maps
    base = fullfile(folder, sprintf('%s_%s', m.method, matlab.lang.makeValidName(m.target)));
    for w = what
        switch w{1}
            case 'stat', v = m.stat;
            case 'p', if all(isnan(m.p)), continue; end, v = m.p;
            case 'fdr', if all(isnan(m.p)), continue; end, v = m.stat .* m.fdr_mask;
            case 'fwer', if isnan(m.fwer_threshold), continue; end, v = m.stat .* m.fwer_mask;
        end
        f = [base '_' w{1} '.nii'];
        save_volume(v, T, f);
        files{end + 1, 1} = f; %#ok<AGROW>
    end
end
end
