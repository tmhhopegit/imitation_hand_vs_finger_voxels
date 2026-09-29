function [z, p] = vlsm_ranksum(lesions, y)
%VLSM_RANKSUM Wilcoxon rank-sum test at every voxel (from WXVLSM / MakeWXMap).
%   [z, p] = vlsm_ranksum(lesions, y)
%   z > 0 means patients without a lesion at the voxel rank higher (the
%   sign of ranksum's zval with the unlesioned group first, as MakeWXMap).
%   Normal approximation with tie and continuity corrections, computed for
%   all voxels at once. MATLAB's ranksum switches to an exact p-value when
%   a group is very small (under 10 patients with fewer than 20 in total),
%   so p-values can differ slightly for such voxels.
%   Voxels where either group is empty give NaN.
y = y(:);
N = numel(y);
r = tiedrank(y);
L = double(lesions);
n1 = sum(L, 1)';
n0 = N - n1;
W0 = sum(r) - L' * r;                        % rank sum of the unlesioned
E = n0 * (N + 1) / 2;
[~, ~, grp] = unique(y);
t = accumarray(grp, 1);
tie = sum(t .^ 3 - t) / (N * (N - 1));
sd = sqrt(n0 .* n1 / 12 .* ((N + 1) - tie));
d = W0 - E;
z = (d - 0.5 * sign(d)) ./ sd;
p = 2 * normcdf(-abs(z));
bad = n1 < 1 | n0 < 1 | sd == 0;
z(bad) = NaN; p(bad) = NaN;
end
