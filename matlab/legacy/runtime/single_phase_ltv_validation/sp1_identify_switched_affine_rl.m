function result = sp1_identify_switched_affine_rl(data, segmentation, opts)
%SP1_IDENTIFY_SWITCHED_AFFINE_RL Identify V0+R+L in ON intervals.
%
% ON model:
%
%   v = V0 + R i + L di/dt.
%
% Integrated form used for numerical stability:
%
%   int(v) = V0 (t-t0) + R int(i) + L (i-i0).

if nargin < 3
    opts = struct();
end

conditionLimit = nv3_get_option(opts, 'conditionLimit', 1e8);
minSamples = nv3_get_option(opts, 'minSamples', 20);

segments = segmentation.onSegments;
rows = cell(height(segments), 13);
details = struct('indices', {}, 'theta', {}, 'fit', {}, 'predictionIntegral', {});
for s = 1:height(segments)
    idx = segments.StartIndex(s):segments.EndIndex(s);
    t = data.t(idx);
    v = data.v(idx);
    i = data.i(idx);
    if numel(idx) < minSamples
        [theta, fit, residual] = invalid_fit();
        passive = false;
        valid = false;
        prediction = NaN(size(t));
    else
        tau = t - t(1);
        vInt = cumtrapz(t, v);
        iInt = cumtrapz(t, i);
        dI = i - i(1);
        X = [tau, iInt, dI];
        [theta, fit] = normalized_solve(X, vInt);
        prediction = X * theta;
        residual = norm(vInt - prediction) / max(norm(vInt), eps);
        passive = all(isfinite(theta)) && theta(2) > 0 && theta(3) > 0;
        valid = passive && fit.rank == 3 && ...
            fit.conditionNumber <= conditionLimit && isfinite(residual);
    end

    rows(s, :) = {s, segments.StartTime(s), segments.EndTime(s), ...
        segments.Duration(s), segments.Samples(s), theta(1), theta(2), ...
        theta(3), fit.rank, fit.conditionNumber, residual, passive, valid};
    details(s).indices = idx;
    details(s).theta = theta;
    details(s).fit = fit;
    details(s).predictionIntegral = prediction;
end

estimates = cell2table(rows, 'VariableNames', { ...
    'Segment', 'StartTime', 'EndTime', 'Duration', 'Samples', ...
    'V0', 'R', 'L', 'Rank', 'ConditionNumber', 'IntegralResidual', ...
    'Passive', 'Valid'});

result = struct();
result.model = 'switched_affine_rl';
result.estimates = estimates;
result.summary = summarize_estimates(estimates, segmentation);
result.details = details;
result.segmentation = segmentation;

end

function [theta, fit, residual] = invalid_fit()
theta = [NaN; NaN; NaN];
fit = struct('rank', 0, 'conditionNumber', Inf, 'singularValues', NaN);
residual = Inf;
end

function [theta, fit] = normalized_solve(X, y)
mask = all(isfinite(X), 2) & isfinite(y);
Xv = X(mask, :);
yv = y(mask);
m = size(X, 2);
theta = NaN(m, 1);
fit = struct('rank', 0, 'conditionNumber', Inf, 'singularValues', NaN);
if size(Xv, 1) < m
    return;
end
colScale = sqrt(mean(Xv .^ 2, 1)).';
yScale = sqrt(mean(yv .^ 2));
colScale(~isfinite(colScale) | colScale <= eps) = 1;
yScale = max(yScale, eps);
Xs = Xv ./ colScale.';
ys = yv ./ yScale;
sv = svd(Xs, 'econ');
fit.singularValues = sv(:).';
if isempty(sv) || max(sv) <= eps
    return;
end
tol = max(size(Xs)) * eps(max(sv));
fit.rank = sum(sv > tol);
if numel(sv) >= m && sv(m) > tol
    fit.conditionNumber = sv(1) / sv(m);
end
if fit.rank < m
    return;
end
gamma = Xs \ ys;
theta = gamma .* yScale ./ colScale;
end

function summary = summarize_estimates(estimates, segmentation)
valid = estimates.Valid > 0;
summary = struct();
summary.SegmentCount = height(estimates);
summary.ValidSegmentCoverage = mean(valid);
summary.OnCoverage = segmentation.stats.OnCoverage;
summary.OffCoverage = segmentation.stats.OffCoverage;
summary.TransitionCoverage = segmentation.stats.TransitionCoverage;
summary.MedianV0 = median_finite(estimates.V0(valid));
summary.MedianR = median_finite(estimates.R(valid));
summary.MedianL = median_finite(estimates.L(valid));
summary.RobustCvV0 = robust_cv(estimates.V0(valid));
summary.RobustCvR = robust_cv(estimates.R(valid));
summary.RobustCvL = robust_cv(estimates.L(valid));
summary.MedianResidual = median_finite(estimates.IntegralResidual(valid));
summary.MedianCondition = median_finite(estimates.ConditionNumber(valid));
end

function value = robust_cv(x)
x = x(isfinite(x));
if isempty(x)
    value = Inf;
else
    value = iqr_local(x) / max(abs(median(x)), eps);
end
end

function value = iqr_local(x)
x = sort(x(:));
value = percentile(x, 75) - percentile(x, 25);
end

function value = percentile(x, p)
if isempty(x)
    value = NaN;
    return;
end
pos = 1 + (numel(x) - 1) * p / 100;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    value = x(lo);
else
    value = x(lo) + (x(hi) - x(lo)) * (pos - lo);
end
end

function value = median_finite(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = median(x);
end
end
