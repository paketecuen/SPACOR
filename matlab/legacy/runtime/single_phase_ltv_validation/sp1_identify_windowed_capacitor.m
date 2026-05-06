function result = sp1_identify_windowed_capacitor(data, config, opts)
%SP1_IDENTIFY_WINDOWED_CAPACITOR Local pure-capacitor equivalent.
%
% In each window:
%
%   i ~= C_w dv/dt.
%
% This handles nearly pure capacitive laboratory cases where adding a
% conductance column makes G numerically arbitrary and the raw power residual
% becomes dominated by measurement noise.

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

rows = cell(numel(starts), 12);
details = struct('indices', {}, 'theta', {}, 'fit', {}, 'signals', {});
for w = 1:numel(starts)
    idx = starts(w):(starts(w) + windowSamples - 1);
    t = data.t(idx);
    vFit = nv3_harmonic_fit(data.v(idx), t, data.f0, harmonics, ...
        struct('maxOrder', 1));
    iFit = nv3_harmonic_fit(data.i(idx), t, data.f0, harmonics, ...
        struct('maxOrder', 0));

    signals = struct('t', t, 'v', vFit.values{1}, ...
        'dv', vFit.values{2}, 'i', iFit.values{1});
    X = signals.dv;
    y = signals.i;
    [theta, fit] = normalized_solve(X, y);
    prediction = X * theta;
    powerResidual = norm(data.v(idx).*data.i(idx) - ...
        signals.v .* prediction) / max(norm(data.v(idx).*data.i(idx)), eps);
    passive = isfinite(theta) && theta > 0;
    valid = passive && fit.rank == 1 && fit.conditionNumber <= conditionLimit && ...
        isfinite(fit.relativeResidual);

    rows(w, :) = {w, mean(t), t(1), t(end), mat2str(harmonics), ...
        theta, fit.rank, fit.conditionNumber, fit.relativeResidual, ...
        powerResidual, passive, valid};

    details(w).indices = idx;
    details(w).theta = theta;
    details(w).fit = fit;
    details(w).signals = signals;
    details(w).prediction = prediction;
end

windows = cell2table(rows, 'VariableNames', { ...
    'Window', 'CenterTime', 'StartTime', 'EndTime', 'Harmonics', ...
    'C', 'Rank', 'ConditionNumber', 'EquationResidual', ...
    'PowerResidual', 'Passive', 'Valid'});

result = struct();
result.model = 'windowed_capacitor_lti';
result.windows = windows;
result.summary = summarize_windows(windows);
result.details = details;
result.config = config;
result.config.windowCycles = windowCycles;
result.config.hopCycles = hopCycles;
result.config.harmonics = harmonics;

end

function [theta, fit] = normalized_solve(X, y)
mask = isfinite(X) & isfinite(y);
Xv = X(mask);
yv = y(mask);
theta = NaN;
fit = struct('rank', 0, 'conditionNumber', Inf, 'relativeResidual', Inf, ...
    'singularValues', NaN);
if numel(Xv) < 1
    return;
end
xScale = sqrt(mean(Xv .^ 2));
yScale = sqrt(mean(yv .^ 2));
xScale = max(xScale, eps);
yScale = max(yScale, eps);
Xs = Xv ./ xScale;
ys = yv ./ yScale;
sv = norm(Xs);
fit.singularValues = sv;
if sv <= eps
    return;
end
fit.rank = 1;
fit.conditionNumber = 1;
gamma = Xs \ ys;
theta = gamma .* yScale ./ xScale;
fit.relativeResidual = norm(yv - Xv * theta) / max(norm(yv), eps);
end

function summary = summarize_windows(windows)
valid = windows.Valid > 0;
summary = struct();
summary.WindowCount = height(windows);
summary.ValidCoverage = mean(valid);
summary.MedianC = median_finite(windows.C(valid));
summary.MedianCondition = median_finite(windows.ConditionNumber(valid));
summary.MedianEquationResidual = median_finite(windows.EquationResidual(valid));
summary.MedianPowerResidual = median_finite(windows.PowerResidual(valid));
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
