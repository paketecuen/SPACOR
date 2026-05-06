function result = sp1_identify_windowed_parallel_rlc(data, config, opts)
%SP1_IDENTIFY_WINDOWED_PARALLEL_RLC Local LTI parallel R-L-C approximation.
%
% In each window:
%
%   i ~= G_w v + Gamma_w primitive(v) + C_w dv/dt,   Gamma_w = 1/L_w.

if nargin < 2 || isempty(config)
    config = sp1_default_config();
end
if nargin < 3
    opts = struct();
end

windowCycles = nv3_get_option(opts, 'windowCycles', config.windowCycles);
hopCycles = nv3_get_option(opts, 'hopCycles', config.hopCycles);
harmonics = nv3_get_option(opts, 'harmonics', config.harmonics);
conditionLimit = nv3_get_option(opts, 'conditionLimit', config.conditionLimit);

windowSamples = max(12, round(windowCycles * data.fs / data.f0));
hopSamples = max(1, round(hopCycles * data.fs / data.f0));
n = numel(data.t);
starts = 1:hopSamples:(n - windowSamples + 1);

rows = cell(numel(starts), 16);
details = struct('indices', {}, 'theta', {}, 'fit', {}, 'signals', {});
for w = 1:numel(starts)
    idx = starts(w):(starts(w) + windowSamples - 1);
    t = data.t(idx);
    vFit = nv3_harmonic_fit(data.v(idx), t, data.f0, harmonics, ...
        struct('maxOrder', 1));
    iFit = nv3_harmonic_fit(data.i(idx), t, data.f0, harmonics, ...
        struct('maxOrder', 0));

    signals = struct('t', t, 'v', vFit.values{1}, ...
        'qv', vFit.primitive, 'dv', vFit.values{2}, ...
        'i', iFit.values{1});
    X = [signals.v, signals.qv, signals.dv];
    y = signals.i;
    [theta, fit] = normalized_solve(X, y);
    prediction = X * theta;
    powerResidual = norm(data.v(idx).*data.i(idx) - ...
        signals.v .* prediction) / max(norm(data.v(idx).*data.i(idx)), eps);
    passive = all(isfinite(theta)) && theta(1) > 0 && theta(2) > 0 && ...
        theta(3) > 0;
    valid = passive && fit.rank == 3 && fit.conditionNumber <= conditionLimit && ...
        isfinite(fit.relativeResidual) && isfinite(powerResidual);
    R = 1 / theta(1);
    L = 1 / theta(2);

    rows(w, :) = {w, mean(t), t(1), t(end), mat2str(harmonics), ...
        theta(1), R, theta(2), L, theta(3), fit.rank, ...
        fit.conditionNumber, fit.relativeResidual, powerResidual, ...
        passive, valid};

    details(w).indices = idx;
    details(w).theta = theta;
    details(w).fit = fit;
    details(w).signals = signals;
    details(w).prediction = prediction;
end

windows = cell2table(rows, 'VariableNames', { ...
    'Window', 'CenterTime', 'StartTime', 'EndTime', 'Harmonics', ...
    'G', 'R', 'Gamma', 'L', 'C', 'Rank', 'ConditionNumber', ...
    'EquationResidual', 'PowerResidual', 'Passive', 'Valid'});

result = struct();
result.model = 'windowed_parallel_rlc_lti';
result.windows = windows;
result.summary = summarize_windows(windows);
result.details = details;
result.config = config;
result.config.windowCycles = windowCycles;
result.config.hopCycles = hopCycles;
result.config.harmonics = harmonics;

end

function [theta, fit] = normalized_solve(X, y)
mask = all(isfinite(X), 2) & isfinite(y);
Xv = X(mask, :);
yv = y(mask);
m = size(X, 2);
theta = NaN(m, 1);
fit = struct('rank', 0, 'conditionNumber', Inf, 'relativeResidual', Inf, ...
    'singularValues', NaN);
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
fit.relativeResidual = norm(yv - Xv * theta) / max(norm(yv), eps);
end

function summary = summarize_windows(windows)
valid = windows.Valid > 0;
summary = struct();
summary.WindowCount = height(windows);
summary.ValidCoverage = mean(valid);
summary.MedianG = median_finite(windows.G(valid));
summary.MedianR = median_finite(windows.R(valid));
summary.MedianGamma = median_finite(windows.Gamma(valid));
summary.MedianL = median_finite(windows.L(valid));
summary.MedianC = median_finite(windows.C(valid));
summary.MedianCondition = median_finite(windows.ConditionNumber(valid));
summary.MedianEquationResidual = median_finite(windows.EquationResidual(valid));
summary.MedianPowerResidual = median_finite(windows.PowerResidual(valid));
summary.RobustCvG = sp1_robust_cv(windows.G(valid));
summary.RobustCvR = sp1_robust_cv(windows.R(valid));
summary.RobustCvGamma = sp1_robust_cv(windows.Gamma(valid));
summary.RobustCvL = sp1_robust_cv(windows.L(valid));
summary.RobustCvC = sp1_robust_cv(windows.C(valid));
end

function value = median_finite(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = median(x);
end
end
