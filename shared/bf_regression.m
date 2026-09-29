function logBF = bf_regression(R2, n, k, prior)
%BF_REGRESSION Log Bayes factor for a linear regression against intercept only.
%   logBF = bf_regression(R2, n, k, prior)
%   R2: R-squared of a model with k predictors (any array); n: observations.
%   logBF > 0 favours an effect; |logBF| > log(3) ~ 1.1 is usually called
%   "substantial", > log(10) ~ 2.3 "strong".
%   prior: 'jzs' (default): Zellner-Siow / JZS prior (Liang et al. 2008;
%            Rouder & Morey 2012), integrated numerically. With k = 1 this
%            is exactly the JZS correlation test of Wetzels & Wagenmakers
%            (2012) that jzs_corbf computed (R2 = r^2).
%          'unit_information': Zellner's g-prior with g = n, in closed form
%            (the prior bayes_anova1 used). bayes_anova1 raised (1 + g) to
%            (n - k - 1)/2 with k = number of GROUPS; the exponent needs the
%            number of predictors (groups - 1), so its log BFs were too low
%            by 0.5*log(n + 1), e.g. 2.3 with 100 patients.
if nargin < 4 || isempty(prior), prior = 'jzs'; end
R2 = min(max(R2, 0), 1 - 1e-12);
switch prior
    case 'unit_information'
        g = n;
        logBF = ((n - 1 - k) / 2) * log1p(g) - ((n - 1) / 2) * log1p(g * (1 - R2));
    case 'jzs'
        % integrate over u = log(g) on a fixed grid (the integrand is smooth
        % in u), in chunks to limit memory; log-sum-exp for stability
        u = linspace(-15, 30, 3001)';
        du = u(2) - u(1);
        g = exp(u);
        w = ones(size(u)); w([1 end]) = 0.5;             % trapezoid weights
        logprior = 0.5 * log(n / 2) - gammaln(0.5) - 0.5 * u - n ./ (2 * g);   % log(pi(g) * g)
        base = ((n - 1 - k) / 2) * log1p(g) + logprior + log(w * du);
        logBF = nan(size(R2));
        chunk = 5000;
        for s = 1:chunk:numel(R2)
            idx = s:min(s + chunk - 1, numel(R2));
            L = base - ((n - 1) / 2) * log1p(g * (1 - R2(idx)'));
            m = max(L, [], 1);
            logBF(idx) = m + log(sum(exp(L - m), 1));
        end
    otherwise
        error('imitation:prior', 'prior must be ''jzs'' or ''unit_information''.');
end
end
