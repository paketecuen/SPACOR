function result = sp1_identify_windowed_rl(data, config, opts)
%SP1_IDENTIFY_WINDOWED_RL Local LTI approximation for R + saturating L.
%
% In each window:
%   v ~= R_w i + L_w di/dt
%
% This is intentionally a local constant-parameter approximation to an LTV
% device.  It should be judged by residual and by whether L_w follows the
% underlying incremental inductance trend, not by exact constant recovery.

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
referenceR = nv3_get_option(opts, 'referenceR', truth_r(data));

windowSamples = max(12, round(windowCycles * data.fs / data.f0));
hopSamples = max(1, round(hopCycles * data.fs / data.f0));
n = numel(data.t);
starts = 1:hopSamples:(n - windowSamples + 1);

rows = cell(numel(starts), 26);
details = struct('indices', {}, 'theta', {}, 'fit', {}, 'signals', {});
for w = 1:numel(starts)
    idx = starts(w):(starts(w) + windowSamples - 1);
    t = data.t(idx);
    vFit = nv3_harmonic_fit(data.v(idx), t, data.f0, harmonics, ...
        struct('maxOrder', 0));
    iFit = nv3_harmonic_fit(data.i(idx), t, data.f0, harmonics, ...
        struct('maxOrder', 1));

    signals = struct('t', t, 'v', vFit.values{1}, ...
        'i', iFit.values{1}, 'di', iFit.values{2});
    X = [signals.i, signals.di];
    y = signals.v;
    [theta, fit] = normalized_solve(X, y);
    prediction = X * theta;
    powerResidual = norm(data.v(idx).*data.i(idx) - prediction.*signals.i) / ...
        max(norm(data.v(idx).*data.i(idx)), eps);
    passive = all(isfinite(theta)) && theta(1) > 0 && theta(2) > 0;
    valid = passive && fit.rank == 2 && fit.conditionNumber <= conditionLimit && ...
        isfinite(fit.relativeResidual) && isfinite(powerResidual);

    truthLTrace = truth_l_trace(data);
    truthDldtTrace = truth_dldt_trace(data);
    truthLMean = mean(truthLTrace(idx));
    truthLMedian = median(truthLTrace(idx));
    truthDldtMean = mean(truthDldtTrace(idx), 'omitnan');
    if isfinite(referenceR) && isfinite(truthDldtMean)
        truthREffMean = referenceR + truthDldtMean;
    else
        truthREffMean = NaN;
    end
    if isfinite(referenceR)
        rLoss = theta(1) - referenceR;
    else
        rLoss = NaN;
    end
    [truthTheta, truthFit] = truth_projection(data, idx);
    truthProjectedR = truthTheta(1);
    truthProjectedL = truthTheta(2);
    if isfinite(referenceR)
        truthProjectedRLoss = truthProjectedR - referenceR;
    else
        truthProjectedRLoss = NaN;
    end
    truthAbsIMean = truth_abs_i_mean(data, idx);
    relLToMean = abs(theta(2) - truthLMean) / max(abs(truthLMean), eps);
    relRlossToMean = abs(rLoss - truthDldtMean) / ...
        max(abs(truthDldtMean), eps);
    relLToProjected = abs(theta(2) - truthProjectedL) / ...
        max(abs(truthProjectedL), eps);
    relRlossToProjected = abs(rLoss - truthProjectedRLoss) / ...
        max(abs(truthProjectedRLoss), eps);

    rows(w, :) = {w, mean(t), t(1), t(end), mat2str(harmonics), ...
        theta(1), theta(2), rLoss, truthLMean, truthLMedian, ...
        truthDldtMean, truthREffMean, truthProjectedR, truthProjectedL, ...
        truthProjectedRLoss, truthAbsIMean, relLToMean, relRlossToMean, ...
        relLToProjected, relRlossToProjected, ...
        fit.rank, fit.conditionNumber, fit.relativeResidual, powerResidual, ...
        passive, valid};

    details(w).indices = idx;
    details(w).theta = theta;
    details(w).fit = fit;
    details(w).signals = signals;
    details(w).prediction = prediction;
    details(w).relativeLToMean = relLToMean;
    details(w).truthProjection = struct('theta', truthTheta, ...
        'fit', truthFit, 'relativeL', relLToProjected, ...
        'relativeRLoss', relRlossToProjected);
end

windows = cell2table(rows, 'VariableNames', { ...
    'Window', 'CenterTime', 'StartTime', 'EndTime', 'Harmonics', ...
    'R', 'L', 'RLoss', 'TruthLMean', 'TruthLMedian', 'TruthDldtMean', ...
    'TruthREffMean', 'TruthProjectedR', 'TruthProjectedL', ...
    'TruthProjectedRLoss', 'TruthAbsIMean', 'RelativeLToMean', ...
    'RelativeRlossToMean', 'RelativeLToProjected', ...
    'RelativeRlossToProjected', 'Rank', 'ConditionNumber', ...
    'EquationResidual', 'PowerResidual', 'Passive', 'Valid'});

result = struct();
result.model = 'windowed_rl_lti';
result.windows = windows;
result.summary = summarize_windows(windows, referenceR);
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

function summary = summarize_windows(windows, referenceR)
valid = windows.Valid > 0;
summary = struct();
summary.WindowCount = height(windows);
summary.ValidCoverage = mean(valid);
summary.MedianR = median_finite(windows.R(valid));
summary.MedianL = median_finite(windows.L(valid));
summary.MedianTruthLMean = median_finite(windows.TruthLMean(valid));
summary.MedianCondition = median_finite(windows.ConditionNumber(valid));
summary.MedianEquationResidual = median_finite(windows.EquationResidual(valid));
summary.MedianPowerResidual = median_finite(windows.PowerResidual(valid));
if isfinite(referenceR)
    summary.RelativeRError = abs(summary.MedianR - referenceR) / ...
        max(abs(referenceR), eps);
else
    summary.RelativeRError = NaN;
end
summary.RelativeLToMedianTruth = abs(summary.MedianL - summary.MedianTruthLMean) / ...
    max(abs(summary.MedianTruthLMean), eps);
summary.MedianRLoss = median_finite(windows.RLoss(valid));
summary.MedianTruthDldtMean = median_finite(windows.TruthDldtMean(valid));
summary.RelativeRLossToTruth = abs(summary.MedianRLoss - ...
    summary.MedianTruthDldtMean) / max(abs(summary.MedianTruthDldtMean), eps);
summary.MedianTruthProjectedR = median_finite(windows.TruthProjectedR(valid));
summary.MedianTruthProjectedL = median_finite(windows.TruthProjectedL(valid));
summary.MedianTruthProjectedRLoss = ...
    median_finite(windows.TruthProjectedRLoss(valid));
summary.RelativeLToProjected = abs(summary.MedianL - ...
    summary.MedianTruthProjectedL) / max(abs(summary.MedianTruthProjectedL), eps);
summary.RelativeRLossToProjected = abs(summary.MedianRLoss - ...
    summary.MedianTruthProjectedRLoss) / ...
    max(abs(summary.MedianTruthProjectedRLoss), eps);
end

function R = truth_r(data)
if isfield(data, 'truth') && isfield(data.truth, 'R')
    R = data.truth.R;
else
    R = NaN;
end
end

function [theta, fit] = truth_projection(data, idx)
theta = [NaN; NaN];
fit = struct('rank', 0, 'conditionNumber', Inf, 'relativeResidual', Inf, ...
    'singularValues', NaN);
if isfield(data, 'clean') && isfield(data.clean, 'i') && ...
        isfield(data.clean, 'di') && isfield(data.clean, 'v')
    X = [data.clean.i(idx), data.clean.di(idx)];
    y = data.clean.v(idx);
    [theta, fit] = normalized_solve(X, y);
end
end

function value = truth_abs_i_mean(data, idx)
if isfield(data, 'clean') && isfield(data.clean, 'i')
    value = mean(abs(data.clean.i(idx)));
else
    value = mean(abs(data.i(idx)));
end
end

function values = truth_l_trace(data)
if isfield(data, 'clean') && isfield(data.clean, 'Lcoef')
    values = data.clean.Lcoef;
elseif isfield(data, 'clean') && isfield(data.clean, 'Linc')
    values = data.clean.Linc;
elseif isfield(data, 'clean') && isfield(data.clean, 'L')
    values = data.clean.L;
else
    values = NaN(size(data.t));
end
end

function values = truth_dldt_trace(data)
if isfield(data, 'clean') && isfield(data.clean, 'dLdt')
    values = data.clean.dLdt;
else
    values = NaN(size(data.t));
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
