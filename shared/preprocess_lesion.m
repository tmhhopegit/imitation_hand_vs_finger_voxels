function X = preprocess_lesion(X, opts)
%PREPROCESS_LESION Optional clean-up of a 3-D binary lesion image (from GetLesions).
%   X = preprocess_lesion(X, opts)
%   opts.gap_fill_slices  s > 0: fill gaps between lesion drawings made on
%                         every few axial slices. A voxel becomes lesioned
%                         if the triangular-weighted count of lesioned
%                         voxels within s-1 slices above and below reaches
%                         s/2 (GetLesions used s = 6). 0 = off.
%   opts.one_hemisphere   true: mirror left-right if most of the lesion is
%                         in the first half of the x axis, then clear that
%                         half (x indices below the midline), so every
%                         lesion lies in the same hemisphere. GetLesions did
%                         this for 91-voxel-wide images (midline x = 46).
if isfield(opts, 'gap_fill_slices') && opts.gap_fill_slices > 0
    s = opts.gap_fill_slices;
    K = reshape([1:s, (s - 1):-1:1], 1, 1, []);
    X = X | convn(double(X), K, 'same') >= s / 2;
end
if isfield(opts, 'one_hemisphere') && opts.one_hemisphere
    mid = ceil(size(X, 1) / 2);                      % 46 for 91 voxels
    if nnz(X(1:mid - 1, :, :)) > nnz(X(mid:end, :, :))
        X = flip(X, 1);
    end
    X(1:mid - 1, :, :) = false;
end
end
