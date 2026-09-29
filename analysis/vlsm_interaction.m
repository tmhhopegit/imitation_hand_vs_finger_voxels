function [p, b] = vlsm_interaction(lesions, y, factors)
%VLSM_INTERACTION Does the effect of a lesion depend on a factor? (MakeInteractionMap_Thresh)
%   [p, b] = vlsm_interaction(lesions, y, factors)
%   At each voxel fits  y ~ lesion + factors + lesion x factors(:,1)  and
%   returns the p-value and coefficient of the interaction term. With one
%   degree of freedom this p equals the type-III ANOVA p that anovan gave
%   for the same model (MakeANOVAmap), but it is computed directly.
%   y: N x 1 (e.g. improvement); factors: N x F, continuous, first column
%   is the one that interacts with the lesion.
%   Voxels where the model cannot be estimated give NaN.
y = y(:);
N = numel(y);
F = [ones(N, 1), factors];
V = size(lesions, 2);
p = nan(V, 1); b = nan(V, 1);
k = size(F, 2) + 2;
for v = 1:V
    l = double(lesions(:, v));
    if all(l == l(1)), continue; end
    X = [F, l, l .* factors(:, 1)];
    [Q, R] = qr(X, 0);
    if abs(R(end, end)) < 1e-10 * abs(R(1, 1)), continue; end
    beta = R \ (Q' * y);
    e = y - X * beta;
    s2 = (e' * e) / (N - k);
    Ri = R \ eye(k);
    se = sqrt(s2 * sum(Ri(end, :) .^ 2));
    b(v) = beta(end);
    p(v) = 2 * tcdf(-abs(beta(end) / se), N - k);
end
end
