function [threshold, significant] = fdr_threshold(p, q, method)
%FDR_THRESHOLD False discovery rate threshold for a set of p-values.
%   [threshold, significant] = fdr_threshold(p, q, method)
%   p: any array; NaNs (voxels not tested) are ignored, NOT counted as tests.
%   method: 'bh' (Benjamini & Hochberg 1995; fdr_bh's 'pdep', default) or
%           'by' (Benjamini & Yekutieli 2001, any dependence; fdr_bh's 'dep').
%   threshold: the largest p-value declared significant (0 if none);
%   significant: p <= threshold (false where p is NaN).
%
%   The Achilles scripts passed whole 91x109x91 p-maps to fdr_bh, whose
%   untested voxels held p = 0. Those ~850,000 "p = 0 tests" counted as
%   both tests and discoveries, which pushed the threshold to almost the
%   uncorrected level (about 0.94*q, or 0.94*q/c(m) for 'dep'). Only tested
%   voxels count here.
if nargin < 3, method = 'bh'; end
v = sort(p(~isnan(p(:))));
m = numel(v);
threshold = 0;
found = false;
if m > 0
    switch method
        case 'bh', c = 1;
        case 'by', c = sum(1 ./ (1:m));
        otherwise, error('imitation:fdr', 'method must be ''bh'' or ''by''.');
    end
    k = find(v <= (1:m)' * q / (m * c), 1, 'last');
    if ~isempty(k), threshold = v(k); found = true; end
end
significant = found & p <= threshold & ~isnan(p);
end
