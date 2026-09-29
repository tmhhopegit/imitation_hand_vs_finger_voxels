function d = vlsm_lpd(lesions, impaired)
%VLSM_LPD Lesion proportion difference at every voxel (Rudrauf et al. 2008).
%   d = vlsm_lpd(lesions, impaired)
%   d = (share of impaired patients lesioned at the voxel)
%     - (share of unimpaired patients lesioned there).
%   Positive values: damage there is more common in impaired patients.
impaired = logical(impaired(:));
L = double(lesions);
d = (L' * impaired) / max(nnz(impaired), 1) - (L' * ~impaired) / max(nnz(~impaired), 1);
end
