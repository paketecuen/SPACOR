function result = sp1_identify_windowed_parallel_series_rl(data, config, opts)
%SP1_IDENTIFY_WINDOWED_PARALLEL_SERIES_RL Windowed R || (Rs+L).
%
% Circuit:
%
%   i = Gp v + i_s
%   v = Rs i_s + L d(i_s)/dt
%
% Eliminating the internal branch current gives
%
%   i = a v + b di/dt + c dv/dt
%
% with
%
%   a = 1/Rs + Gp,   b = -L/Rs,   c = L Gp/Rs.
%
% The default numerical path fits the integrated form
%
%   primitive(i) = a primitive(v) + b i + c v + k
%
% in each window.  This avoids differentiating measured signals directly and
% makes the geometric rank deficiency of pure sinusoidal windows explicit.

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

rows = cell(numel(starts), 19);
details = struct('indices', {}, 'theta', {}, 'fit', {}, 'signals', {});
for w = 1:numel(starts)
    idx = starts(w):(starts(w) + windowSamples - 1);
    t = data.t(idx);
    vFit = nv3_harmonic_fit(data.v(idx), t, data.f0, harmonics, ...
        struct('maxOrder', 1));
    iFit = nv3_harmonic_fit(data.i(idx), t, data.f0, harmonics, ...
        struct('maxOrder', 1));

    signals = struct('t', t, ...
        'v', vFit.values{1}, 'dv', vFit.values{2}, ...
        'qv', vFit.primitive, ...
        'i', iFit.values{1}, 'di', iFit.values{2}, ...
        'qi', iFit.primitive);
    X = [signals.qv, signals.i, signals.v, ones(size(t))];
    y = signals.qi;
    [theta, fit] = normalized_solve(X, y);

    a = theta(1);
    b = theta(2);
    c = theta(3);
    intercept = theta(4);
    [Gp, Rp, Rs, L] = physical_parameters(a, b, c);
    qiPrediction = X * theta;
    iPrediction = a * signals.v + b * signals.di + c * signals.dv;
    powerResidual = norm(data.v(idx).*data.i(idx) - ...
        signals.v .* iPrediction) / max(norm(data.v(idx).*data.i(idx)), eps);
    passive = all(isfinite([Gp, Rp, Rs, L])) && ...
        Gp > 0 && Rp > 0 && Rs > 0 && L > 0;
    valid = passive && fit.rank == 4 && fit.conditionNumber <= conditionLimit && ...
        isfinite(fit.relativeResidual) && isfinite(powerResidual);

    rows(w, :) = {w, mean(t), t(1), t(end), mat2str(harmonics), ...
        a, b, c, intercept, Gp, Rp, Rs, L, fit.rank, ...
        fit.conditionNumber, fit.relativeResidual, powerResidual, ...
        passive, valid};

    details(w).indices = idx;
    details(w).theta = theta;
    details(w).fit = fit;
    details(w).signals = signals;
    details(w).qiPrediction = qiPrediction;
    details(w).iPrediction = iPrediction;
end

windows = cell2table(rows, 'VariableNames', { ...
    'Window', 'CenterTime', 'StartTime', 'EndTime', 'Harmonics', ...
    'A', 'B', 'Ccoef', 'Intercept', 'Gp', 'Rp', 'Rs', 'L', ...
    'Rank', 'ConditionNumber', 'EquationResidual', 'PowerResidual', ...
    'Passive', 'Valid'});

result = struct();
result.model = 'windowed_parallel_series_rl_lti';
result.windows = windows;
result.summary = summarize_windows(windows);
result.details = details;
result.config = config;
result.config.windowCycles = windowCycles;
result.config.hopCycles = hopCycles;
result.config.harmonics = harmonics;

end

function [Gp, Rp, Rs, L] = physical_parameters(a, b, c)
Gp = NaN;
Rp = NaN;
Rs = NaN;
L = NaN;
if ~all(isfinite([a, b, c])) || abs(b) <= eps
    return;
end
Gp = -c / b;
den = a - Gp;
if ~isfinite(Gp) || ~isfinite(den) || abs(den) <= eps
    return;
end
Rs = 1 / den;
L = -b * Rs;
if isfinite(Gp) && abs(Gp) > eps
    Rp = 1 / Gp;
end
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
summary.MedianA = median_finite(windows.A(valid));
summary.MedianB = median_finite(windows.B(valid));
summary.MedianCcoef = median_finite(windows.Ccoef(valid));
summary.MedianGp = median_finite(windows.Gp(valid));
summary.MedianRp = median_finite(windows.Rp(valid));
summary.MedianRs = median_finite(windows.Rs(valid));
summary.MedianL = median_finite(windows.L(valid));
summary.MedianCondition = median_finite(windows.ConditionNumber(valid));
summary.MedianEquationResidual = median_finite(windows.EquationResidual(valid));
summary.MedianPowerResidual = median_finite(windows.PowerResidual(valid));
summary.RobustCvGp = sp1_robust_cv(windows.Gp(valid));
summary.RobustCvRp = sp1_robust_cv(windows.Rp(valid));
summary.RobustCvRs = sp1_robust_cv(windows.Rs(valid));
summary.RobustCvL = sp1_robust_cv(windows.L(valid));
end

function value = median_finite(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = median(x);
end
end
