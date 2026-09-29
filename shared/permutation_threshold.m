function [threshold, null_max] = permutation_threshold(stat_fun, Y, n_perm, alpha, seed, workers)
%PERMUTATION_THRESHOLD Family-wise error threshold by permuting patients.
%   [threshold, null_max] = permutation_threshold(stat_fun, Y, n_perm, alpha, seed, workers)
%   stat_fun(Yp) returns one column of statistics per map for the permuted
%   scores Yp (rows = patients); larger |statistic| = stronger effect.
%   Each permutation shuffles whole ROWS of Y, so every score of a patient
%   moves together and difference maps keep their pairing. BatesVLSM,
%   FrequentistSingleRegionComparison and LPDD permuted the two tasks
%   separately, which breaks the pairing and gives the wrong null for the
%   difference between tasks.
%   The same permutations are used for every call with the same seed.
%   threshold(m): the (1 - alpha) quantile of the maximum |statistic| over
%   voxels, per map; null_max: n_perm x maps.
if nargin < 4 || isempty(alpha), alpha = 0.05; end
if nargin < 5 || isempty(seed), seed = 1; end
if nargin < 6 || isempty(workers), workers = 0; end
N = size(Y, 1);
stream = RandStream('mt19937ar', 'Seed', seed);
orders = zeros(N, n_perm);
for k = 1:n_perm, orders(:, k) = randperm(stream, N); end
m = size(stat_fun(Y), 2);
null_max = zeros(n_perm, m);
parfor (k = 1:n_perm, workers)
    S = stat_fun(Y(orders(:, k), :));
    null_max(k, :) = max(abs(S), [], 1, 'omitnan');
end
sorted = sort(null_max, 1, 'descend');
threshold = sorted(max(1, ceil(n_perm * alpha)), :);
end
