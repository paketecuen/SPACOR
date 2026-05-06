function result = sp1_identify_switched_lti_library(data, segmentation, opts)
%SP1_IDENTIFY_SWITCHED_LTI_LIBRARY Identify canonical switched ON models.
%
% The switch is treated as a visibility/topology state.  In ON intervals we
% test canonical series-port equivalents:
%
%   switched_R_like:   v = V0 + R i
%   switched_RL_like:  v = V0 + R i + L di/dt
%   switched_RC_like:  v = V0 + R i + Gamma int(i)dt
%   switched_RLC_like: v = V0 + R i + L di/dt + Gamma int(i)dt
%
% The fit is performed in integral form so the ON-state library does not need
% numerical differentiation:
%
%   int(v) = V0 tau + R int(i) + L (i-i0) + Gamma int(int(i)).

if nargin < 3
    opts = struct();
end

models = nv3_get_option(opts, 'models', ...
    {'switched_R_like', 'switched_RL_like', ...
    'switched_RC_like', 'switched_RLC_like'});
conditionLimit = nv3_get_option(opts, 'conditionLimit', 1e8);
minSamples = nv3_get_option(opts, 'minSamples', 20);

segments = segmentation.onSegments;
rows = {};
details = struct('Model', {}, 'Segment', {}, 'indices', {}, ...
    'theta', {}, 'fit', {}, 'predictionIntegral', {});
detailIdx = 0;
for modelIdx = 1:numel(models)
    modelName = char(models{modelIdx});
    spec = model_spec(modelName);
    for s = 1:height(segments)
        idx = segments.StartIndex(s):segments.EndIndex(s);
        t = data.t(idx);
        v = data.v(idx);
        i = data.i(idx);
        if numel(idx) < minSamples
            [thetaFull, fit, residual] = invalid_fit();
            passive = false;
            valid = false;
            prediction = NaN(size(t));
        else
            [X, y] = integral_system(spec, t, v, i);
            [theta, fit] = normalized_solve(X, y);
            prediction = X * theta;
            residual = norm(y - prediction) / max(norm(y), eps);
            thetaFull = expand_theta(spec, theta);
            passive = is_passive(spec, thetaFull);
            valid = passive && fit.rank == numel(spec.terms) && ...
                fit.conditionNumber <= conditionLimit && isfinite(residual);
        end

        rows(end + 1, :) = {modelName, spec.complexity, s, ...
            segments.StartTime(s), segments.EndTime(s), ...
            segments.Duration(s), segments.Samples(s), thetaFull(1), ...
            thetaFull(2), thetaFull(3), thetaFull(4), ...
            safe_inverse(thetaFull(4)), fit.rank, fit.conditionNumber, ...
            residual, passive, valid}; %#ok<AGROW>
        detailIdx = detailIdx + 1;
        details(detailIdx).Model = modelName;
        details(detailIdx).Segment = s;
        details(detailIdx).indices = idx;
        details(detailIdx).theta = thetaFull;
        details(detailIdx).fit = fit;
        details(detailIdx).predictionIntegral = prediction;
    end
end

if isempty(rows)
    estimatesAll = empty_estimates();
else
    estimatesAll = cell2table(rows, 'VariableNames', variable_names());
end

modelTable = summarize_models(estimatesAll, segmentation, opts);
if isempty(modelTable) || height(modelTable) == 0
    best = unresolved_summary(segmentation);
    bestEstimates = empty_estimates();
else
    best = select_best_model(modelTable, opts);
    bestEstimates = estimatesAll(strcmp(estimatesAll.Model, best.Model), :);
end

result = struct();
result.model = 'switching_lti_library';
result.bestModel = best.Model;
result.estimates = bestEstimates;
result.allEstimates = estimatesAll;
result.modelTable = modelTable;
result.summary = best;
result.details = details;
result.segmentation = segmentation;

end

function best = unresolved_summary(segmentation)
best = struct();
best.Model = 'switched_unresolved';
best.PhysicalMeaning = ['switching-like activity detected, but ON segments ' ...
    'are not parameter-identifiable'];
best.Accepted = false;
best.Confidence = 0;
best.OnCoverage = segmentation.stats.OnCoverage;
best.OffCoverage = segmentation.stats.OffCoverage;
best.TransitionCoverage = segmentation.stats.TransitionCoverage;
best.ValidSegmentCoverage = 0;
best.MedianResidual = Inf;
best.StabilityCV = Inf;
best.Complexity = 0;
best.MedianV0 = NaN;
best.MedianR = NaN;
best.MedianL = NaN;
best.MedianGamma = NaN;
best.MedianC = NaN;
best.RobustCvR = Inf;
best.RobustCvL = Inf;
best.RobustCvGamma = Inf;
best.MedianCondition = Inf;
best.Reason = sprintf(['on %.2f, off %.2f, transition %.2f; no valid ' ...
    'ON core segments for parameter identification'], ...
    segmentation.stats.OnCoverage, segmentation.stats.OffCoverage, ...
    segmentation.stats.TransitionCoverage);
end

function spec = model_spec(name)
switch name
    case 'switched_R'
        name = 'switched_R_like';
        terms = {'V0', 'R'};
        physical = 'switched R-like ON equivalent: v = V0 + R i';
    case 'switched_R_like'
        terms = {'V0', 'R'};
        physical = 'switched R-like ON equivalent: v = V0 + R i';
    case 'switched_RL'
        name = 'switched_RL_like';
        terms = {'V0', 'R', 'L'};
        physical = 'switched R-L-like ON equivalent: v = V0 + R i + L di/dt';
    case 'switched_RL_like'
        terms = {'V0', 'R', 'L'};
        physical = 'switched R-L-like ON equivalent: v = V0 + R i + L di/dt';
    case 'switched_RC'
        name = 'switched_RC_like';
        terms = {'V0', 'R', 'Gamma'};
        physical = ['switched R-C-like ON equivalent: v = V0 + R i + ' ...
            'Gamma int(i)dt'];
    case 'switched_RC_like'
        terms = {'V0', 'R', 'Gamma'};
        physical = ['switched R-C-like ON equivalent: v = V0 + R i + ' ...
            'Gamma int(i)dt'];
    case 'switched_RLC'
        name = 'switched_RLC_like';
        terms = {'V0', 'R', 'L', 'Gamma'};
        physical = ['switched R-L-C-like ON equivalent: v = V0 + R i + ' ...
            'L di/dt + Gamma int(i)dt'];
    case 'switched_RLC_like'
        terms = {'V0', 'R', 'L', 'Gamma'};
        physical = ['switched R-L-C-like ON equivalent: v = V0 + R i + L di/dt + ' ...
            'Gamma int(i)dt'];
    otherwise
        error('Unknown switched ON model "%s".', name);
end
spec = struct('name', name, 'terms', {terms}, ...
    'complexity', numel(terms), 'physicalMeaning', physical);
end

function [X, y] = integral_system(spec, t, v, i)
t = t(:);
v = v(:);
i = i(:);
tau = t - t(1);
vInt = cumtrapz(t, v);
iInt = cumtrapz(t, i);
deltaI = i - i(1);
iiInt = cumtrapz(t, iInt);
y = vInt;
X = zeros(numel(t), numel(spec.terms));
for k = 1:numel(spec.terms)
    switch spec.terms{k}
        case 'V0'
            X(:, k) = tau;
        case 'R'
            X(:, k) = iInt;
        case 'L'
            X(:, k) = deltaI;
        case 'Gamma'
            X(:, k) = iiInt;
        otherwise
            error('Unknown switched term "%s".', spec.terms{k});
    end
end
end

function thetaFull = expand_theta(spec, theta)
thetaFull = [NaN; NaN; NaN; NaN]; % V0, R, L, Gamma
for k = 1:numel(spec.terms)
    switch spec.terms{k}
        case 'V0'
            thetaFull(1) = theta(k);
        case 'R'
            thetaFull(2) = theta(k);
        case 'L'
            thetaFull(3) = theta(k);
        case 'Gamma'
            thetaFull(4) = theta(k);
    end
end
end

function passive = is_passive(spec, thetaFull)
passive = all(isfinite(thetaFull([1, 2])) | isnan(thetaFull([1, 2])));
for k = 1:numel(spec.terms)
    switch spec.terms{k}
        case 'R'
            passive = passive && isfinite(thetaFull(2)) && thetaFull(2) > 0;
        case 'L'
            passive = passive && isfinite(thetaFull(3)) && thetaFull(3) > 0;
        case 'Gamma'
            passive = passive && isfinite(thetaFull(4)) && thetaFull(4) > 0;
    end
end
end

function summaries = summarize_models(estimates, segmentation, opts)
modelValues = string(estimates.Model);
models = unique(modelValues, 'stable');
rows = cell(numel(models), 20);
for k = 1:numel(models)
    name = char(models(k));
    W = estimates(modelValues == models(k), :);
    valid = W.Valid > 0;
    spec = model_spec(name);
    stability = model_stability(spec, W, valid);
    coverageOk = segmentation.stats.OnCoverage >= ...
        nv3_get_option(opts, 'onCoverageGate', 0.10) && ...
        segmentation.stats.OffCoverage >= ...
        nv3_get_option(opts, 'offCoverageGate', 0.10);
    validCoverage = mean(valid);
    validOk = validCoverage >= ...
        nv3_get_option(opts, 'validSegmentCoverageGate', 0.80);
    residual = median_finite(W.IntegralResidual(valid));
    residualOk = residual <= ...
        nv3_get_option(opts, 'integralResidualGate', 2e-2);
    stableOk = stability <= nv3_get_option(opts, 'stabilityGate', 0.20);
    accepted = coverageOk && validOk && residualOk && stableOk;
    confidence = mean([double(coverageOk), double(validOk), ...
        double(residualOk), double(stableOk)]);
    reason = sprintf(['on %.2f, off %.2f, valid segments %.2f, ', ...
        'residual %.3g, stabilityCV %.3g'], ...
        segmentation.stats.OnCoverage, segmentation.stats.OffCoverage, ...
        validCoverage, residual, stability);

    rows(k, :) = {string(name), string(spec.physicalMeaning), ...
        accepted, confidence, ...
        segmentation.stats.OnCoverage, segmentation.stats.OffCoverage, ...
        segmentation.stats.TransitionCoverage, validCoverage, ...
        residual, stability, spec.complexity, ...
        median_finite(W.V0(valid)), median_finite(W.R(valid)), ...
        median_finite(W.L(valid)), median_finite(W.Gamma(valid)), ...
        median_finite(W.C(valid)), ...
        robust_cv(W.R(valid)), robust_cv(W.L(valid)), ...
        robust_cv(W.Gamma(valid)), reason};
end

summaries = cell2table(rows, 'VariableNames', { ...
    'Model', 'PhysicalMeaning', 'Accepted', 'Confidence', 'OnCoverage', ...
    'OffCoverage', 'TransitionCoverage', 'ValidSegmentCoverage', ...
    'MedianResidual', 'StabilityCV', 'Complexity', 'MedianV0', ...
    'MedianR', 'MedianL', 'MedianGamma', 'MedianC', 'RobustCvR', ...
    'RobustCvL', 'RobustCvGamma', 'Reason'});
end

function best = select_best_model(modelTable, opts)
accepted = find(modelTable.Accepted > 0);
if isempty(accepted)
    [~, idx] = max(modelTable.Confidence);
else
    residuals = modelTable.MedianResidual(accepted);
    [bestResidual, localBest] = min(residuals);
    bestIdx = accepted(localBest);
    factor = nv3_get_option(opts, 'switchingEquivalenceFactor', 1.5);
    absTol = nv3_get_option(opts, 'switchingEquivalenceAbsTol', 2e-3);
    near = accepted(modelTable.MedianResidual(accepted) <= ...
        max(bestResidual * factor, bestResidual + absTol));
    [~, simpleLocal] = min(modelTable.Complexity(near));
    idx = near(simpleLocal);
    if modelTable.Complexity(bestIdx) < modelTable.Complexity(idx)
        idx = bestIdx;
    end
end

best = struct();
best.Model = char(modelTable.Model(idx));
best.PhysicalMeaning = char(modelTable.PhysicalMeaning(idx));
best.Accepted = modelTable.Accepted(idx);
best.Confidence = modelTable.Confidence(idx);
best.OnCoverage = modelTable.OnCoverage(idx);
best.OffCoverage = modelTable.OffCoverage(idx);
best.TransitionCoverage = modelTable.TransitionCoverage(idx);
best.ValidSegmentCoverage = modelTable.ValidSegmentCoverage(idx);
best.MedianResidual = modelTable.MedianResidual(idx);
best.StabilityCV = modelTable.StabilityCV(idx);
best.Complexity = modelTable.Complexity(idx);
best.MedianV0 = modelTable.MedianV0(idx);
best.MedianR = modelTable.MedianR(idx);
best.MedianL = modelTable.MedianL(idx);
best.MedianGamma = modelTable.MedianGamma(idx);
best.MedianC = modelTable.MedianC(idx);
best.RobustCvR = modelTable.RobustCvR(idx);
best.RobustCvL = modelTable.RobustCvL(idx);
best.RobustCvGamma = modelTable.RobustCvGamma(idx);
best.Reason = char(modelTable.Reason(idx));
best.MedianCondition = NaN;
end

function stability = model_stability(spec, W, valid)
values = [];
for k = 1:numel(spec.terms)
    switch spec.terms{k}
        case 'R'
            values(end + 1) = robust_cv(W.R(valid)); %#ok<AGROW>
        case 'L'
            values(end + 1) = robust_cv(W.L(valid)); %#ok<AGROW>
        case 'Gamma'
            values(end + 1) = robust_cv(W.Gamma(valid)); %#ok<AGROW>
    end
end
if isempty(values)
    stability = Inf;
else
    stability = max(values);
end
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

function estimates = empty_estimates()
estimates = cell2table(cell(0, numel(variable_names())), ...
    'VariableNames', variable_names());
end

function names = variable_names()
names = {'Model', 'Complexity', 'Segment', 'StartTime', 'EndTime', ...
    'Duration', 'Samples', 'V0', 'R', 'L', 'Gamma', 'C', 'Rank', ...
    'ConditionNumber', 'IntegralResidual', 'Passive', 'Valid'};
end

function [theta, fit, residual] = invalid_fit()
theta = [NaN; NaN; NaN; NaN];
fit = struct('rank', 0, 'conditionNumber', Inf, 'singularValues', NaN);
residual = Inf;
end

function value = safe_inverse(x)
if isfinite(x) && x > 0
    value = 1 / x;
else
    value = NaN;
end
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
