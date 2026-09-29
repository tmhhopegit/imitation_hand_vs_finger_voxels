function logBF = vlsm_bayes(lesions, y, prior)
%VLSM_BAYES Log Bayes factor for a lesioned-vs-not difference at every voxel.
%   logBF = vlsm_bayes(lesions, y)            JZS prior (default)
%   logBF = vlsm_bayes(lesions, y, prior)     'jzs' or 'unit_information'
%   logBF > 0 favours a difference; < 0 favours no difference (log(1/3) =
%   -1.1 is "substantial" evidence for none). The Bayesian maps were meant
%   to show where the data support NO effect, which p-values cannot.
%
%   This is a one-predictor regression of y on the lesion indicator
%   (equivalent to the one-way ANOVA in bayes_anova1), computed for all
%   voxels at once. bayes_anova1 used the unit-information prior with an
%   exponent that was off by one (see bf_regression): its maps were lower
%   everywhere by 0.5*log(n+1). Use prior 'unit_information' to get the
%   corrected version of the original.
%   Voxels where everyone (or no one) is lesioned give NaN.
if nargin < 3, prior = 'jzs'; end
y = y(:) - mean(y);
N = numel(y);
L = double(lesions);
n1 = sum(L, 1)';
n0 = N - n1;
s1 = L' * y;                                  % sum over lesioned (overall mean is 0)
SST = sum(y .^ 2);
% between-group sum of squares with the grand mean at 0: s1^2/n1 + s0^2/n0, s0 = -s1
R2 = (s1 .^ 2 ./ n1 + s1 .^ 2 ./ n0) / SST;
logBF = bf_regression(R2, N, 1, prior);
logBF(n1 < 1 | n0 < 1) = NaN;
end
