function result = nv3_identify_three_wire_windowed(data, config, opts)
%NV3_IDENTIFY_THREE_WIRE_WINDOWED Windowed matrix identification for 3-wire.
%
% Implemented models:
%   delta_parallel_gl : delta branches with parallel G and Gamma=1/L.
%   wye_series_rl     : wye branches with series R and L.

if nargin < 2 || isempty(config)
    config = nv3_default_config();
end
if nargin < 3
    opts = struct();
end

model = char(nv3_get_option(opts, 'model', data.model));
fs = nv3_get_option(opts, 'fs', data.fs);
f0 = nv3_get_option(opts, 'f0', data.f0);
windowCycles = nv3_get_option(opts, 'windowCycles', config.windowCycles);
hopCycles = nv3_get_option(opts, 'hopCycles', config.hopCycles);
windowSamples = max(12, round(windowCycles * fs / f0));
hopSamples = max(1, round(hopCycles * fs / f0));
conditionLimit = nv3_get_option(opts, 'conditionLimit', config.conditionLimit);

n = numel(data.t);
starts = 1:hopSamples:(n - windowSamples + 1);
paramNames = model_parameter_names(model);

rows = cell(numel(starts), 12 + numel(paramNames));
details = struct('indices', {}, 'theta', {}, 'thetaNames', {}, 'harmonics', {}, ...
    'fit', {}, 'prediction', {});

for w = 1:numel(starts)
    idx = starts(w):(starts(w) + windowSamples - 1);
    t = data.t(idx);
    [signals, selected, fitBundle] = fit_window_signals(data, idx, t, f0, config, opts, model);
    [X, y, thetaNames] = build_model_system(model, signals);
    [theta, fit] = solve_normalized_ls(X, y);
    params = theta_to_params(model, theta);
    passive = params_passive(params, paramNames);
    prediction = predict_window(model, params, signals);
    powerResidual = power_residual(data, idx, prediction);
    valid = fit.rank >= numel(theta) && fit.conditionNumber <= conditionLimit && ...
        passive && isfinite(fit.relativeResidual) && isfinite(powerResidual);

    rowValues = [{w, mean(t), t(1), t(end), model, mat2str(selected), ...
        fit.rank, fit.conditionNumber, fit.relativeResidual, powerResidual, ...
        passive, valid}, params_to_cells(params, paramNames)];
    rows(w, :) = rowValues;

    details(w).indices = idx;
    details(w).theta = theta;
    details(w).thetaNames = thetaNames;
    details(w).harmonics = selected;
    details(w).fit = fit;
    details(w).fitBundle = fitBundle;
    details(w).prediction = prediction;
end

varNames = [{'Window', 'CenterTime', 'StartTime', 'EndTime', 'Model', 'Harmonics', ...
    'Rank', 'ConditionNumber', 'EquationResidual', 'PowerResidual', ...
    'Passive', 'Valid'}, paramNames];
windows = cell2table(rows, 'VariableNames', varNames);

result = struct();
result.model = model;
result.windows = windows;
result.details = details;
result.summary = summarize_windows(windows, data.truth, paramNames);
result.config = config;

end

function names = model_parameter_names(model)
def = nv3_model_definition(model);
switch def.topology
    case 'delta_parallel'
        names = def.parameterNames;
    case 'wye_series'
        if strcmp(model, 'wye_series_rl')
            names = {'Ra', 'Rb', 'Rc', 'La', 'Lb', 'Lc'};
        else
            names = def.parameterNames;
        end
    otherwise
        error('Unknown three-wire model "%s".', model);
end
end

function [signals, selected, fitBundle] = fit_window_signals(data, idx, t, f0, config, opts, model)
raw = [data.vab(idx), data.vbc(idx), data.vca(idx), ...
    data.ia(idx), data.ib(idx), data.ic(idx)];

harmonics = nv3_get_option(opts, 'harmonics', config.harmonics);
adaptive = nv3_get_option(opts, 'adaptiveHarmonics', config.adaptiveHarmonics);
if ischar(harmonics) || isstring(harmonics)
    adaptive = strcmpi(char(harmonics), 'adaptive');
end
if adaptive
    candidates = nv3_get_option(opts, 'harmonicCandidates', config.harmonicCandidates);
    selection = nv3_select_active_harmonics(raw, t, f0, candidates, merge_options(config, opts));
    selected = selection.selectedHarmonics;
else
    selected = harmonics;
end

def = nv3_model_definition(model);
voltageOrder = double(any(strcmp(def.terms, 'C')) && strcmp(def.topology, 'delta_parallel'));
fitOpts0 = struct('maxOrder', voltageOrder);
fitOpts1 = struct('maxOrder', 1);
vabFit = nv3_harmonic_fit(data.vab(idx), t, f0, selected, fitOpts0);
vbcFit = nv3_harmonic_fit(data.vbc(idx), t, f0, selected, fitOpts0);
vcaFit = nv3_harmonic_fit(data.vca(idx), t, f0, selected, fitOpts0);
iaFit = nv3_harmonic_fit(data.ia(idx), t, f0, selected, fitOpts1);
ibFit = nv3_harmonic_fit(data.ib(idx), t, f0, selected, fitOpts1);
icFit = nv3_harmonic_fit(data.ic(idx), t, f0, selected, fitOpts1);

if strcmp(def.topology, 'delta_parallel')
    % Enforce KVL in line voltages and their primitive/derivative coordinates.
    vca = -vabFit.values{1} - vbcFit.values{1};
    % Enforce KVL in the primitive coordinates to avoid numerical drift.
    qca = -vabFit.primitive - vbcFit.primitive;
    if voltageOrder >= 1
        dvca = -vabFit.values{2} - vbcFit.values{2};
    else
        dvca = [];
    end
else
    vca = vcaFit.values{1};
    qca = vcaFit.primitive;
    dvca = [];
end

signals = struct();
signals.t = t;
signals.vab = vabFit.values{1};
signals.vbc = vbcFit.values{1};
signals.vca = vca;
signals.qab = vabFit.primitive;
signals.qbc = vbcFit.primitive;
signals.qca = qca;
if voltageOrder >= 1
    signals.dvab = vabFit.values{2};
    signals.dvbc = vbcFit.values{2};
    signals.dvca = dvca;
end
signals.ia = iaFit.values{1};
signals.ib = ibFit.values{1};
signals.ic = icFit.values{1};
signals.dia = iaFit.values{2};
signals.dib = ibFit.values{2};
signals.dic = icFit.values{2};

fitBundle = struct('vab', vabFit, 'vbc', vbcFit, 'vca', vcaFit, ...
    'ia', iaFit, 'ib', ibFit, 'ic', icFit);
end

function optsOut = merge_options(config, opts)
optsOut = opts;
fields = {'activeThreshold', 'maxBasisCondition'};
for idx = 1:numel(fields)
    name = fields{idx};
    if ~isfield(optsOut, name)
        optsOut.(name) = config.(name);
    end
end
end

function [X, y, names] = build_model_system(model, s)
def = nv3_model_definition(model);
switch def.topology
    case 'delta_parallel'
        X = [];
        for termIdx = 1:numel(def.terms)
            feature = delta_term_features(def.terms{termIdx}, s);
            z = zeros(size(feature.ab));
            Xterm = [ ...
                feature.ab, z, -feature.ca; ...
                -feature.ab, feature.bc, z; ...
                z, -feature.bc, feature.ca];
            X = [X, Xterm]; %#ok<AGROW>
        end
        y = [s.ia; s.ib; s.ic];
        names = def.thetaNames;

    case 'wye_series'
        if ~strcmp(model, 'wye_series_rl')
            error('Windowed wye-series model "%s" is not implemented yet.', model);
        end
        X = [ ...
            s.ia, -s.ib, zeros(size(s.ia)), s.dia, -s.dib, zeros(size(s.ia)); ...
            zeros(size(s.ia)), s.ib, -s.ic, zeros(size(s.ia)), s.dib, -s.dic; ...
            -s.ia, zeros(size(s.ia)), s.ic, -s.dia, zeros(size(s.ia)), s.dic];
        y = [s.vab; s.vbc; s.vca];
        names = {'Ra', 'Rb', 'Rc', 'La', 'Lb', 'Lc'};
end
end

function feature = delta_term_features(term, s)
switch term
    case 'G'
        feature = struct('ab', s.vab, 'bc', s.vbc, 'ca', s.vca);
    case 'Gamma'
        feature = struct('ab', s.qab, 'bc', s.qbc, 'ca', s.qca);
    case 'C'
        feature = struct('ab', s.dvab, 'bc', s.dvbc, 'ca', s.dvca);
    otherwise
        error('Unknown delta term "%s".', term);
end
end

function [theta, fit] = solve_normalized_ls(X, y)
mask = all(isfinite(X), 2) & isfinite(y);
Xv = X(mask, :);
yv = y(mask);
m = size(X, 2);
theta = NaN(m, 1);
fit = struct('rank', 0, 'conditionNumber', Inf, 'relativeResidual', Inf, ...
    'validRows', sum(mask), 'rowCount', numel(mask));
if size(Xv, 1) < m
    return;
end
colScale = sqrt(mean(Xv .^ 2, 1)).';
yScale = sqrt(mean(yv .^ 2));
colScale(~isfinite(colScale) | colScale <= eps) = 1;
yScale = max(yScale, eps);
Xs = Xv ./ colScale.';
ys = yv ./ yScale;
s = svd(Xs, 'econ');
tol = max(size(Xs)) * eps(max(s));
fit.rank = sum(s > tol);
if numel(s) >= m && s(m) > tol
    fit.conditionNumber = s(1) / s(m);
end
if fit.rank < m
    return;
end
gamma = Xs \ ys;
theta = gamma .* yScale ./ colScale;
fit.relativeResidual = norm(yv - Xv * theta) / max(norm(yv), eps);
end

function params = theta_to_params(model, theta)
def = nv3_model_definition(model);
switch def.topology
    case 'delta_parallel'
        params = struct();
        for idx = 1:numel(def.thetaNames)
            params.(def.thetaNames{idx}) = theta(idx);
        end
        if isfield(params, 'Gab')
            params.Rab = 1 / params.Gab;
            params.Rbc = 1 / params.Gbc;
            params.Rca = 1 / params.Gca;
        end
        if isfield(params, 'Gammaab')
            params.Lab = 1 / params.Gammaab;
            params.Lbc = 1 / params.Gammabc;
            params.Lca = 1 / params.Gammaca;
        end
    case 'wye_series'
        if ~strcmp(model, 'wye_series_rl')
            error('Windowed wye-series model "%s" is not implemented yet.', model);
        end
        params = struct('Ra', theta(1), 'Rb', theta(2), 'Rc', theta(3), ...
            'La', theta(4), 'Lb', theta(5), 'Lc', theta(6));
end
end

function cells = params_to_cells(params, names)
cells = cell(1, numel(names));
for idx = 1:numel(names)
    name = names{idx};
    if isfield(params, name)
        cells{idx} = params.(name);
    else
        cells{idx} = NaN;
    end
end
end

function ok = params_passive(params, names)
ok = true;
for idx = 1:numel(names)
    name = names{idx};
    value = params.(name);
    ok = ok && isfinite(value) && value > 0;
end
end

function prediction = predict_window(model, params, s)
def = nv3_model_definition(model);
switch def.topology
    case 'delta_parallel'
        iab = zeros(size(s.vab));
        ibc = zeros(size(s.vbc));
        ica = zeros(size(s.vca));
        for termIdx = 1:numel(def.terms)
            term = def.terms{termIdx};
            feature = delta_term_features(term, s);
            iab = iab + params.([term 'ab']) * feature.ab;
            ibc = ibc + params.([term 'bc']) * feature.bc;
            ica = ica + params.([term 'ca']) * feature.ca;
        end
        ia = iab - ica;
        ib = ibc - iab;
        ic = ica - ibc;
        prediction = struct('vab', s.vab, 'vbc', s.vbc, 'vca', s.vca, ...
            'ia', ia, 'ib', ib, 'ic', ic, 'iab', iab, 'ibc', ibc, 'ica', ica);
    case 'wye_series'
        if ~strcmp(model, 'wye_series_rl')
            error('Windowed wye-series model "%s" is not implemented yet.', model);
        end
        van = params.Ra * s.ia + params.La * s.dia;
        vbn = params.Rb * s.ib + params.Lb * s.dib;
        vcn = params.Rc * s.ic + params.Lc * s.dic;
        prediction = struct('vab', van - vbn, 'vbc', vbn - vcn, ...
            'vca', vcn - van, 'ia', s.ia, 'ib', s.ib, 'ic', s.ic, ...
            'van', van, 'vbn', vbn, 'vcn', vcn);
end
prediction.power = prediction.vab(:) .* prediction.ia(:) + ...
    prediction.vbc(:) .* (prediction.ia(:) + prediction.ib(:));
end

function residual = power_residual(data, idx, prediction)
pMeasured = data.vab(idx) .* data.ia(idx) + data.vbc(idx) .* ...
    (data.ia(idx) + data.ib(idx));
pPredicted = prediction.power;
residual = norm(pMeasured(:) - pPredicted(:)) / max(norm(pMeasured(:)), eps);
end

function summary = summarize_windows(windows, truth, paramNames)
valid = windows.Valid > 0;
parameter = paramNames(:);
truthValue = NaN(numel(paramNames), 1);
medianEstimate = NaN(numel(paramNames), 1);
relativeError = NaN(numel(paramNames), 1);
iqrEstimate = NaN(numel(paramNames), 1);
for idx = 1:numel(paramNames)
    name = paramNames{idx};
    values = windows.(name);
    mask = valid & isfinite(values);
    if any(mask)
        medianEstimate(idx) = median(values(mask));
        iqrEstimate(idx) = percentile(values(mask), 75) - percentile(values(mask), 25);
    end
    if isfield(truth, name)
        truthValue(idx) = truth.(name);
        relativeError(idx) = abs(medianEstimate(idx) - truthValue(idx)) / ...
            max(abs(truthValue(idx)), eps);
    end
end
summary = table(parameter, truthValue, medianEstimate, relativeError, iqrEstimate, ...
    'VariableNames', {'Parameter', 'Truth', 'Median', 'RelativeError', 'IQR'});
summary.Properties.UserData.ValidCoverage = mean(valid);
summary.Properties.UserData.MedianCondition = median_finite(windows.ConditionNumber(valid));
summary.Properties.UserData.MedianEquationResidual = median_finite(windows.EquationResidual(valid));
summary.Properties.UserData.MedianPowerResidual = median_finite(windows.PowerResidual(valid));
end

function value = percentile(x, p)
x = sort(x(isfinite(x)));
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
