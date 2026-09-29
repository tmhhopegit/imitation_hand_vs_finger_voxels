function [t, p, df] = vlsm_ttest(lesions, y, vartype)
%VLSM_TTEST Two-sample t-test at every voxel (Bates et al. 2003 VLSM).
%   [t, p, df] = vlsm_ttest(lesions, y)            Welch (unequal variances)
%   [t, p, df] = vlsm_ttest(lesions, y, 'equal')   pooled variance
%   lesions: N x V logical; y: N x 1 scores (no NaNs).
%   t > 0 means patients WITHOUT a lesion at the voxel score higher, i.e.
%   damage there goes with lower scores (the sign MakeTMap used).
%   Voxels with fewer than 2 patients in either group give NaN.
%
%   All voxels are computed at once from group sums, instead of calling
%   ttest2 once per voxel as MakeTMap did (about 100-1000x faster).
if nargin < 3, vartype = 'unequal'; end
y = y(:) - mean(y);                 % centring improves the precision of the sums
N = numel(y);
L = double(lesions);
n1 = sum(L, 1)';                    % lesioned
n0 = N - n1;
s1 = L' * y;          s0 = sum(y) - s1;
q1 = L' * (y .^ 2);   q0 = sum(y .^ 2) - q1;
m1 = s1 ./ n1;        m0 = s0 ./ n0;
v1 = max(q1 - n1 .* m1 .^ 2, 0) ./ (n1 - 1);
v0 = max(q0 - n0 .* m0 .^ 2, 0) ./ (n0 - 1);
if strcmp(vartype, 'equal')
    sp = ((n1 - 1) .* v1 + (n0 - 1) .* v0) ./ (N - 2);
    se = sqrt(sp .* (1 ./ n1 + 1 ./ n0));
    df = repmat(N - 2, size(n1));
else
    a = v0 ./ n0; b = v1 ./ n1;
    se = sqrt(a + b);
    df = (a + b) .^ 2 ./ (a .^ 2 ./ (n0 - 1) + b .^ 2 ./ (n1 - 1));
end
t = (m0 - m1) ./ se;
p = 2 * tcdf(-abs(t), df);
bad = n1 < 2 | n0 < 2 | se == 0;
t(bad) = NaN; p(bad) = NaN; df(bad) = NaN;
end
