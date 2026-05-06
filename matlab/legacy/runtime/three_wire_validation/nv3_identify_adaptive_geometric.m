function result = nv3_identify_adaptive_geometric(data, config, opts)
%NV3_IDENTIFY_ADAPTIVE_GEOMETRIC Adaptive geometric model order for 3-wire RLC.
%
% The default physical family is the three-wire delta parallel RLC model:
%
%   i_xy = G_xy v_xy + Gamma_xy q_xy + C_xy dv_xy/dt
%
% with q_xy = integral(v_xy)dt and Gamma_xy = 1/L_xy.  The routine evaluates a
% hierarchy of symmetry embeddings
%
%   balanced        : 3 free parameters
%   pair_equal      : 6 free parameters
%   full_unbalanced : 9 free parameters
%
% and selects the simplest passive, geometrically adequate model whose
% residuals are close to the best usable candidate in the same window.  In
% black-box mode the candidate library can also include three-wire wye series
% R/L equivalents, so topology is part of the selected model.

if nargin < 2 || isempty(config)
    config = nv3_default_config();
end
if nargin < 3
    opts = struct();
end

config = fill_defaults(config);
model = char(nv3_get_option(opts, 'model', 'delta_parallel_glc'));
if ~strcmp(model, 'delta_parallel_glc')
    error('nv3_identify_adaptive_geometric currently uses delta_parallel_glc coordinates.');
end
candidateMode = char(nv3_get_option(opts, 'candidateMode', 'symmetry_glc'));

windowCycles = nv3_get_option(opts, 'windowCycles', config.windowCycles);
hopCycles = nv3_get_option(opts, 'hopCycles', config.hopCycles);
windowSamples = max(12, round(windowCycles * data.fs / data.f0));
hopSamples = max(1, round(hopCycles * data.fs / data.f0));
n = numel(data.t);
starts = 1:hopSamples:(n - windowSamples + 1);

candidates = candidate_models(candidateMode);
residualGate = nv3_get_option(opts, 'residualGate', infer_residual_gate(data, config));
energyGate = nv3_get_option(opts, 'energyGate', max(config.energyResidualFloor, ...
    config.energyResidualFactor * residualGate));
selectionTolerance = nv3_get_option(opts, 'selectionRelativeTolerance', ...
    config.selectionRelativeTolerance);

windowRows = cell(numel(starts), 23);
candidateTables = cell(numel(starts), 1);
selectedFullTheta = NaN(numel(starts), 9);
details = struct('indices', {}, 'signals', {}, 'evaluations', {}, 'chosen', {});

for w = 1:numel(starts)
    idx = starts(w):(starts(w) + windowSamples - 1);
    signals = fitted_window_signals(data, idx, config, opts);
    evaluations = evaluate_candidates(candidates, signals, data, idx, config);
    chosen = choose_candidate(evaluations, residualGate, energyGate, selectionTolerance);

    if chosen.selected && numel(chosen.fullTheta) == 9
        selectedFullTheta(w, :) = chosen.fullTheta(:).';
    end

    windowRows(w, :) = {w, mean(data.t(idx)), data.t(idx(1)), data.t(idx(end)), ...
        mat2str(signals.harmonics), chosen.name, chosen.topology, chosen.termSet, ...
        chosen.family, chosen.nParameters, chosen.rank, chosen.sigmaMin, ...
        chosen.geometricVolume, chosen.conditionNumber, ...
        chosen.equationResidual, chosen.powerResidual, chosen.energyResidual, ...
        residualGate, energyGate, chosen.passive, chosen.adequate, ...
        chosen.selected, chosen.reason};

    candidateTables{w} = candidate_table(w, mean(data.t(idx)), evaluations, ...
        residualGate, energyGate);

    details(w).indices = idx;
    details(w).signals = signals;
    details(w).evaluations = evaluations;
    details(w).chosen = chosen;
end

windows = cell2table(windowRows, 'VariableNames', { ...
    'Window', 'CenterTime', 'StartTime', 'EndTime', 'Harmonics', ...
    'SelectedModel', 'SelectedTopology', 'SelectedTerms', 'SelectedFamily', ...
    'SelectedN', 'Rank', 'SigmaMin', ...
    'GeometricVolume', 'ConditionNumber', 'EquationResidual', ...
    'PowerResidual', 'EnergyResidual', 'ResidualGate', 'EnergyGate', ...
    'Passive', 'Adequate', 'Selected', 'Reason'});

candidatesTable = vertcat(candidateTables{:});
parameterSummary = summarize_parameters(selectedFullTheta, windows, data);
summary = summarize_windows(windows, selectedFullTheta, data);

result = struct();
result.model = model;
result.windows = windows;
result.candidates = candidatesTable;
result.selectedFullTheta = selectedFullTheta;
result.parameterSummary = parameterSummary;
result.summary = summary;
result.details = details;
result.config = config;
result.residualGate = residualGate;
result.energyGate = energyGate;
result.parameterNames = full_parameter_names();

end

function config = fill_defaults(config)
if ~isfield(config, 'modelResidualFloor')
    config.modelResidualFloor = 2.5e-3;
end
if ~isfield(config, 'modelResidualNoiseFactor')
    config.modelResidualNoiseFactor = 6;
end
if ~isfield(config, 'energyResidualFloor')
    config.energyResidualFloor = 2.5e-3;
end
if ~isfield(config, 'energyResidualFactor')
    config.energyResidualFactor = 5;
end
if ~isfield(config, 'selectionRelativeTolerance')
    config.selectionRelativeTolerance = 0.10;
end
if ~isfield(config, 'rankTolerance')
    config.rankTolerance = 1e-8;
end
if ~isfield(config, 'sigmaFloor')
    config.sigmaFloor = 1e-8;
end
end

function gate = infer_residual_gate(data, config)
noiseGate = 0;
if isfield(data, 'noiseSnrDb') && isfinite(data.noiseSnrDb)
    noiseGate = config.modelResidualNoiseFactor * 10 ^ (-data.noiseSnrDb / 20);
end
gate = max(config.modelResidualFloor, noiseGate);
end

function signals = fitted_window_signals(data, idx, config, opts)
t = data.t(idx);
raw = [data.vab(idx), data.vbc(idx), data.vca(idx), ...
    data.ia(idx), data.ib(idx), data.ic(idx)];

harmonics = nv3_get_option(opts, 'harmonics', config.harmonics);
adaptive = nv3_get_option(opts, 'adaptiveHarmonics', config.adaptiveHarmonics);
if ischar(harmonics) || isstring(harmonics)
    adaptive = strcmpi(char(harmonics), 'adaptive');
end
if adaptive
    selection = nv3_select_active_harmonics(raw, t, data.f0, ...
        config.harmonicCandidates, config);
    h = selection.selectedHarmonics;
else
    h = harmonics;
end

vabFit = nv3_harmonic_fit(data.vab(idx), t, data.f0, h, struct('maxOrder', 1));
vbcFit = nv3_harmonic_fit(data.vbc(idx), t, data.f0, h, struct('maxOrder', 1));
iaFit = nv3_harmonic_fit(data.ia(idx), t, data.f0, h, struct('maxOrder', 1));
ibFit = nv3_harmonic_fit(data.ib(idx), t, data.f0, h, struct('maxOrder', 1));
icFit = nv3_harmonic_fit(data.ic(idx), t, data.f0, h, struct('maxOrder', 1));

signals = struct();
signals.t = t;
signals.vab = vabFit.values{1};
signals.vbc = vbcFit.values{1};
signals.vca = -signals.vab - signals.vbc;
signals.qab = vabFit.primitive;
signals.qbc = vbcFit.primitive;
signals.qca = -signals.qab - signals.qbc;
signals.dvab = vabFit.values{2};
signals.dvbc = vbcFit.values{2};
signals.dvca = -signals.dvab - signals.dvbc;
signals.ia = iaFit.values{1};
signals.ib = ibFit.values{1};
signals.ic = icFit.values{1};
signals.dia = iaFit.values{2};
signals.dib = ibFit.values{2};
signals.dic = icFit.values{2};
signals.qia = iaFit.primitive;
signals.qib = ibFit.primitive;
signals.qic = -signals.qia - signals.qib;
signals.harmonics = h(:).';
signals.fit = struct('vab', vabFit, 'vbc', vbcFit, ...
    'ia', iaFit, 'ib', ibFit, 'ic', icFit);
end

function [X, y] = full_delta_system(s)
z = zeros(size(s.vab));
X = [
    s.vab, z, -s.vca, s.qab, z, -s.qca, s.dvab, z, -s.dvca
    -s.vab, s.vbc, z, -s.qab, s.qbc, z, -s.dvab, s.dvbc, z
    z, -s.vbc, s.vca, z, -s.qbc, s.qca, z, -s.dvbc, s.dvca
    ];
y = [s.ia; s.ib; s.ic];
end

function [X, y] = full_wye_system(s)
z = zeros(size(s.ia));
X = [
    s.ia, -s.ib, z, s.dia, -s.dib, z, s.qia, -s.qib, z
    z, s.ib, -s.ic, z, s.dib, -s.dic, z, s.qib, -s.qic
    -s.ia, z, s.ic, -s.dia, z, s.dic, -s.qia, z, s.qic
    ];
y = [s.vab; s.vbc; s.vca];
end

function candidates = candidate_models(mode)
candidates = struct('name', {}, 'topology', {}, 'termSet', {}, ...
    'family', {}, 'S', {}, 'nParameters', {});
switch mode
    case 'symmetry_glc'
        candidates(end + 1) = make_candidate('balanced', ...
            'delta_parallel', 'glc', ...
            'balanced', symmetry_matrix('balanced'));
        candidates(end + 1) = make_candidate('pair_ab_bc', ...
            'delta_parallel', 'glc', ...
            'pair_equal', symmetry_matrix('pair_ab_bc'));
        candidates(end + 1) = make_candidate('pair_bc_ca', ...
            'delta_parallel', 'glc', ...
            'pair_equal', symmetry_matrix('pair_bc_ca'));
        candidates(end + 1) = make_candidate('pair_ca_ab', ...
            'delta_parallel', 'glc', ...
            'pair_equal', symmetry_matrix('pair_ca_ab'));
        candidates(end + 1) = make_candidate('full_unbalanced', ...
            'delta_parallel', 'glc', ...
            'full_unbalanced', eye(9));
    case {'delta_blackbox', 'blackbox_delta'}
        candidates = delta_blackbox_candidates(candidates);
    case {'mixed_blackbox', 'blackbox_mixed'}
        candidates = delta_blackbox_candidates(candidates);
        candidates = wye_blackbox_candidates(candidates);
    otherwise
        error('Unknown candidate mode "%s".', mode);
end
end

function candidates = delta_blackbox_candidates(candidates)
termSets = {'g', 'gamma', 'c', 'gl', 'gc', 'gammac', 'glc'};
symmetries = {'balanced', 'pair_ab_bc', 'pair_bc_ca', ...
    'pair_ca_ab', 'full_unbalanced'};
for tIdx = 1:numel(termSets)
    termSet = termSets{tIdx};
    for sIdx = 1:numel(symmetries)
        symmetry = symmetries{sIdx};
        S = delta_term_symmetry_matrix(termSet, symmetry);
        name = ['delta_' termSet '_' symmetry];
        family = symmetry_family(symmetry);
        candidates(end + 1) = make_candidate(name, 'delta_parallel', ...
            termSet, family, S); %#ok<AGROW>
    end
end
end

function candidates = wye_blackbox_candidates(candidates)
termSets = {'r', 'l', 'gamma', 'rl', 'rgamma', 'lgamma', 'rlgamma'};
symmetries = {'balanced', 'pair_ab_bc', 'pair_bc_ca', ...
    'pair_ca_ab', 'full_unbalanced'};
for tIdx = 1:numel(termSets)
    termSet = termSets{tIdx};
    for sIdx = 1:numel(symmetries)
        symmetry = symmetries{sIdx};
        S = wye_term_symmetry_matrix(termSet, symmetry);
        name = ['wye_' termSet '_' symmetry];
        family = symmetry_family(symmetry);
        candidates(end + 1) = make_candidate(name, 'wye_series', ...
            termSet, family, S); %#ok<AGROW>
    end
end
end

function candidate = make_candidate(name, topology, termSet, family, S)
candidate = struct('name', name, 'topology', topology, ...
    'termSet', termSet, 'family', family, 'S', S, ...
    'nParameters', size(S, 2));
end

function family = symmetry_family(symmetry)
switch symmetry
    case 'balanced'
        family = 'balanced';
    case {'pair_ab_bc', 'pair_bc_ca', 'pair_ca_ab'}
        family = 'pair_equal';
    case 'full_unbalanced'
        family = 'full_unbalanced';
    otherwise
        error('Unknown symmetry "%s".', symmetry);
end
end

function S = delta_term_symmetry_matrix(termSet, symmetry)
if strcmp(symmetry, 'full_unbalanced')
    base = eye(9);
else
    base = symmetry_matrix(symmetry);
end
termMask = delta_term_indices(termSet);
cols = [];
for col = 1:size(base, 2)
    activeRows = find(base(:, col) ~= 0);
    if any(termMask(activeRows))
        kept = zeros(9, 1);
        kept(activeRows) = base(activeRows, col) .* termMask(activeRows);
        if any(kept ~= 0)
            cols = [cols, kept]; %#ok<AGROW>
        end
    end
end
S = cols;
end

function mask = delta_term_indices(termSet)
mask = false(9, 1);
switch termSet
    case 'g'
        mask(1:3) = true;
    case 'gamma'
        mask(4:6) = true;
    case 'c'
        mask(7:9) = true;
    case 'gl'
        mask(1:6) = true;
    case 'gc'
        mask([1:3, 7:9]) = true;
    case 'gammac'
        mask(4:9) = true;
    case 'glc'
        mask(:) = true;
    otherwise
        error('Unknown term set "%s".', termSet);
end
end

function S = wye_term_symmetry_matrix(termSet, symmetry)
if strcmp(symmetry, 'full_unbalanced')
    base = eye(9);
else
    base = symmetry_matrix(symmetry);
end
termMask = wye_term_indices(termSet);
cols = [];
for col = 1:size(base, 2)
    activeRows = find(base(:, col) ~= 0);
    if any(termMask(activeRows))
        kept = zeros(9, 1);
        kept(activeRows) = base(activeRows, col) .* termMask(activeRows);
        if any(kept ~= 0)
            cols = [cols, kept]; %#ok<AGROW>
        end
    end
end
S = cols;
end

function mask = wye_term_indices(termSet)
mask = false(9, 1);
switch termSet
    case 'r'
        mask(1:3) = true;
    case 'l'
        mask(4:6) = true;
    case 'gamma'
        mask(7:9) = true;
    case 'rl'
        mask(1:6) = true;
    case 'rgamma'
        mask([1:3, 7:9]) = true;
    case 'lgamma'
        mask(4:9) = true;
    case 'rlgamma'
        mask(:) = true;
    otherwise
        error('Unknown wye term set "%s".', termSet);
end
end

function S = symmetry_matrix(kind)
switch kind
    case 'balanced'
        S = [
            1 0 0
            1 0 0
            1 0 0
            0 1 0
            0 1 0
            0 1 0
            0 0 1
            0 0 1
            0 0 1];
    case 'pair_ab_bc'
        S = [
            1 0 0 0 0 0
            1 0 0 0 0 0
            0 1 0 0 0 0
            0 0 1 0 0 0
            0 0 1 0 0 0
            0 0 0 1 0 0
            0 0 0 0 1 0
            0 0 0 0 1 0
            0 0 0 0 0 1];
    case 'pair_bc_ca'
        S = [
            1 0 0 0 0 0
            0 1 0 0 0 0
            0 1 0 0 0 0
            0 0 1 0 0 0
            0 0 0 1 0 0
            0 0 0 1 0 0
            0 0 0 0 1 0
            0 0 0 0 0 1
            0 0 0 0 0 1];
    case 'pair_ca_ab'
        S = [
            1 0 0 0 0 0
            0 1 0 0 0 0
            1 0 0 0 0 0
            0 0 1 0 0 0
            0 0 0 1 0 0
            0 0 1 0 0 0
            0 0 0 0 1 0
            0 0 0 0 0 1
            0 0 0 0 1 0];
    otherwise
        error('Unknown symmetry matrix "%s".', kind);
end
end

function evaluations = evaluate_candidates(candidates, signals, data, idx, config)
evaluations = candidates;
measuredPower = data.vab(idx) .* data.ia(idx) + data.vbc(idx) .* ...
    (data.ia(idx) + data.ib(idx));
for k = 1:numel(candidates)
    switch candidates(k).topology
        case 'delta_parallel'
            [Xfull, y] = full_delta_system(signals);
        case 'wye_series'
            [Xfull, y] = full_wye_system(signals);
        otherwise
            error('Unknown candidate topology "%s".', candidates(k).topology);
    end
    X = Xfull * candidates(k).S;
    [theta, fit] = normalized_solve(X, y, config);
    fullTheta = candidates(k).S * theta;
    passive = all(isfinite(theta)) && all(theta > 0);
    rankOk = fit.rank == candidates(k).nParameters;
    sigmaOk = isfinite(fit.sigmaMin) && fit.sigmaMin > config.sigmaFloor;
    conditionOk = isfinite(fit.conditionNumber) && ...
        fit.conditionNumber <= config.conditionLimit;
    powerResidual = Inf;
    energyResidual = Inf;
    if rankOk && sigmaOk && conditionOk && passive
        switch candidates(k).topology
            case 'delta_parallel'
                prediction = predict_delta(signals, fullTheta);
                energyResidual = delta_energy_residual(signals, measuredPower, ...
                    fullTheta);
            case 'wye_series'
                prediction = predict_wye(signals, fullTheta);
                energyResidual = wye_energy_residual(signals, measuredPower, ...
                    fullTheta);
        end
        powerResidual = norm(measuredPower(:) - prediction.power(:)) / ...
            max(norm(measuredPower(:)), eps);
    end
    evaluations(k).theta = theta;
    evaluations(k).fullTheta = fullTheta;
    evaluations(k).rank = fit.rank;
    evaluations(k).singularValues = fit.singularValues;
    evaluations(k).sigmaMin = fit.sigmaMin;
    evaluations(k).geometricVolume = fit.geometricVolume;
    evaluations(k).conditionNumber = fit.conditionNumber;
    evaluations(k).equationResidual = fit.relativeResidual;
    evaluations(k).powerResidual = powerResidual;
    evaluations(k).energyResidual = energyResidual;
    evaluations(k).passive = passive;
    evaluations(k).geometricValid = rankOk && sigmaOk && conditionOk;
    evaluations(k).usable = rankOk && sigmaOk && conditionOk && passive && ...
        isfinite(powerResidual) && isfinite(energyResidual);
end
end

function chosen = choose_candidate(evaluations, residualGate, energyGate, relativeTolerance)
chosen = empty_choice('reject_no_geometric_candidate');
usableAll = evaluations([evaluations.usable]);
if isempty(usableAll)
    return;
end
allScores = candidate_scores(usableAll, energyGate);
bestScore = min(allScores);
bestOverall = usableAll(find(allScores == bestScore, 1, 'first'));
nValues = unique([usableAll.nParameters]);
for nIdx = 1:numel(nValues)
    group = usableAll([usableAll.nParameters] == nValues(nIdx));
    scores = candidate_scores(group, energyGate);
    [candidateScore, bestIdx] = min(scores);
    candidate = group(bestIdx);
    scoreClose = candidateScore <= (1 + relativeTolerance) * bestScore + ...
        0.10 * residualGate;
    residualOk = candidate.powerResidual <= residualGate && ...
        candidate.equationResidual <= 2 * residualGate && ...
        candidate.energyResidual <= energyGate;
    if residualOk && scoreClose
        chosen = candidate_to_choice(candidate, true, 'adequate_smallest_model');
        return;
    end
end
chosen = candidate_to_choice(bestOverall, false, 'reject_best_candidate_not_adequate');
chosen.name = 'reject';
chosen.topology = 'none';
chosen.termSet = 'none';
chosen.family = 'reject';
chosen.nParameters = 0;
chosen.selected = false;
chosen.adequate = false;
chosen.fullTheta = NaN(9, 1);
end

function scores = candidate_scores(candidates, energyGate)
scores = [candidates.equationResidual] + [candidates.powerResidual] + ...
    min([candidates.energyResidual], energyGate) / max(energyGate, eps);
end

function choice = empty_choice(reason)
choice = struct('name', 'reject', 'family', 'reject', 'nParameters', 0, ...
    'topology', 'none', 'termSet', 'none', ...
    'rank', NaN, 'sigmaMin', NaN, 'geometricVolume', NaN, ...
    'conditionNumber', NaN, 'equationResidual', NaN, ...
    'powerResidual', NaN, 'energyResidual', NaN, 'adequate', false, ...
    'passive', false, 'selected', false, 'reason', reason, ...
    'fullTheta', NaN(9, 1));
end

function choice = candidate_to_choice(candidate, adequate, reason)
choice = struct('name', candidate.name, 'termSet', candidate.termSet, ...
    'topology', candidate.topology, 'family', candidate.family, ...
    'nParameters', candidate.nParameters, 'rank', candidate.rank, ...
    'sigmaMin', candidate.sigmaMin, ...
    'geometricVolume', candidate.geometricVolume, ...
    'conditionNumber', candidate.conditionNumber, ...
    'equationResidual', candidate.equationResidual, ...
    'powerResidual', candidate.powerResidual, ...
    'energyResidual', candidate.energyResidual, ...
    'adequate', adequate, 'passive', candidate.passive, ...
    'selected', adequate, 'reason', reason, 'fullTheta', candidate.fullTheta);
end

function [theta, fit] = normalized_solve(X, y, config)
mask = all(isfinite(X), 2) & isfinite(y);
Xv = X(mask, :);
yv = y(mask);
m = size(X, 2);
theta = NaN(m, 1);
fit = struct('rank', 0, 'conditionNumber', Inf, 'relativeResidual', Inf, ...
    'singularValues', NaN, 'sigmaMin', 0, 'geometricVolume', 0);
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
if isempty(sv) || ~isfinite(max(sv)) || max(sv) <= eps
    return;
end
tol = max(max(size(Xs)) * eps(max(sv)), config.rankTolerance * max(sv));
fit.rank = sum(sv > tol);
if numel(sv) >= m
    fit.sigmaMin = sv(m);
    fit.geometricVolume = prod(sv(1:m));
    if sv(m) > tol
        fit.conditionNumber = sv(1) / sv(m);
    end
end
if fit.rank < m
    return;
end

gamma = Xs \ ys;
theta = gamma .* yScale ./ colScale;
fit.relativeResidual = norm(yv - Xv * theta) / max(norm(yv), eps);
end

function prediction = predict_delta(s, theta)
iab = theta(1) * s.vab + theta(4) * s.qab + theta(7) * s.dvab;
ibc = theta(2) * s.vbc + theta(5) * s.qbc + theta(8) * s.dvbc;
ica = theta(3) * s.vca + theta(6) * s.qca + theta(9) * s.dvca;
ia = iab - ica;
ib = ibc - iab;
ic = ica - ibc;
prediction = struct('iab', iab, 'ibc', ibc, 'ica', ica, ...
    'ia', ia, 'ib', ib, 'ic', ic);
prediction.power = s.vab(:) .* ia(:) + s.vbc(:) .* (ia(:) + ib(:));
end

function prediction = predict_wye(s, theta)
van = theta(1) * s.ia + theta(4) * s.dia;
vbn = theta(2) * s.ib + theta(5) * s.dib;
vcn = theta(3) * s.ic + theta(6) * s.dic;
if numel(theta) >= 9
    van = van + theta(7) * s.qia;
    vbn = vbn + theta(8) * s.qib;
    vcn = vcn + theta(9) * s.qic;
end
vab = van - vbn;
vbc = vbn - vcn;
vca = vcn - van;
prediction = struct('van', van, 'vbn', vbn, 'vcn', vcn, ...
    'vab', vab, 'vbc', vbc, 'vca', vca);
prediction.power = vab(:) .* s.ia(:) + vbc(:) .* (s.ia(:) + s.ib(:));
end

function residual = delta_energy_residual(s, measuredPower, theta)
t = s.t(:);
eIn = trapz(t, measuredPower(:));
pDiss = theta(1) * s.vab(:) .^ 2 + theta(2) * s.vbc(:) .^ 2 + ...
    theta(3) * s.vca(:) .^ 2;
eDiss = trapz(t, pDiss);
stored = 0.5 * theta(4) * s.qab(:) .^ 2 + ...
    0.5 * theta(5) * s.qbc(:) .^ 2 + ...
    0.5 * theta(6) * s.qca(:) .^ 2 + ...
    0.5 * theta(7) * s.vab(:) .^ 2 + ...
    0.5 * theta(8) * s.vbc(:) .^ 2 + ...
    0.5 * theta(9) * s.vca(:) .^ 2;
deltaStored = stored(end) - stored(1);
residual = abs(eIn - eDiss - deltaStored) / max(abs(eIn), eps);
end

function residual = wye_energy_residual(s, measuredPower, theta)
t = s.t(:);
eIn = trapz(t, measuredPower(:));
pDiss = theta(1) * s.ia(:) .^ 2 + theta(2) * s.ib(:) .^ 2 + ...
    theta(3) * s.ic(:) .^ 2;
eDiss = trapz(t, pDiss);
stored = 0.5 * theta(4) * s.ia(:) .^ 2 + ...
    0.5 * theta(5) * s.ib(:) .^ 2 + ...
    0.5 * theta(6) * s.ic(:) .^ 2;
if numel(theta) >= 9
    stored = stored + 0.5 * theta(7) * s.qia(:) .^ 2 + ...
        0.5 * theta(8) * s.qib(:) .^ 2 + ...
        0.5 * theta(9) * s.qic(:) .^ 2;
end
deltaStored = stored(end) - stored(1);
residual = abs(eIn - eDiss - deltaStored) / max(abs(eIn), eps);
end

function tableOut = candidate_table(window, centerTime, evaluations, residualGate, energyGate)
rows = cell(numel(evaluations), 20);
for k = 1:numel(evaluations)
    e = evaluations(k);
    rows(k, :) = {window, centerTime, e.name, e.topology, e.termSet, ...
        e.family, e.nParameters, e.rank, e.sigmaMin, ...
        e.geometricVolume, e.conditionNumber, e.equationResidual, ...
        e.powerResidual, e.energyResidual, residualGate, energyGate, ...
        e.passive, e.geometricValid, e.usable, mat2str(e.fullTheta(:).')};
end
tableOut = cell2table(rows, 'VariableNames', { ...
    'Window', 'CenterTime', 'CandidateModel', 'CandidateTopology', ...
    'CandidateTerms', 'CandidateFamily', 'NParameters', 'Rank', ...
    'SigmaMin', 'GeometricVolume', 'ConditionNumber', ...
    'EquationResidual', 'PowerResidual', 'EnergyResidual', ...
    'ResidualGate', 'EnergyGate', 'Passive', 'GeometricValid', ...
    'Usable', 'FullTheta'});
end

function parameterSummary = summarize_parameters(thetaValues, windows, data)
names = full_parameter_names();
truth = truth_vector(data);
selected = windows.Selected > 0;
rows = cell(numel(names), 6);
for p = 1:numel(names)
    values = thetaValues(selected, p);
    values = values(isfinite(values));
    med = NaN;
    iqrValue = NaN;
    rel = NaN;
    if ~isempty(values)
        med = median(values);
        iqrValue = percentile(values, 75) - percentile(values, 25);
    end
    if isfinite(truth(p)) && isfinite(med)
        rel = abs(med - truth(p)) / max(abs(truth(p)), eps);
    end
    rows(p, :) = {names{p}, truth(p), med, rel, iqrValue, numel(values)};
end
parameterSummary = cell2table(rows, 'VariableNames', { ...
    'Parameter', 'Truth', 'MedianEstimate', 'RelativeError', 'IQR', ...
    'ValidCount'});
end

function summary = summarize_windows(windows, thetaValues, data)
selected = windows.Selected > 0;
adequate = windows.Adequate > 0;
truth = truth_vector(data);
paramMedian = median_theta(thetaValues(selected, :));
rel = abs(paramMedian(:) - truth(:)) ./ max(abs(truth(:)), eps);
summary = struct();
summary.WindowCount = height(windows);
summary.SelectedCoverage = mean(selected);
summary.AdequateCoverage = mean(adequate);
summary.MedianSelectedN = median_finite(windows.SelectedN(selected));
summary.DominantFamily = dominant_family(windows.SelectedFamily(selected));
summary.MedianCondition = median_finite(windows.ConditionNumber(selected));
summary.MedianSigmaMin = median_finite(windows.SigmaMin(selected));
summary.MedianGeometricVolume = median_finite(windows.GeometricVolume(selected));
summary.MedianEquationResidual = median_finite(windows.EquationResidual(selected));
summary.MedianPowerResidual = median_finite(windows.PowerResidual(selected));
summary.MedianEnergyResidual = median_finite(windows.EnergyResidual(selected));
summary.MaxBranchRelativeError = max(rel, [], 'omitnan');
end

function family = dominant_family(values)
if isempty(values)
    family = 'reject';
    return;
end
families = {'balanced', 'pair_equal', 'full_unbalanced', 'reject'};
counts = zeros(size(families));
for k = 1:numel(families)
    counts(k) = sum(strcmp(values, families{k}));
end
[~, idx] = max(counts);
family = families{idx};
end

function theta = truth_vector(data)
names = full_parameter_names();
theta = NaN(numel(names), 1);
if ~isfield(data, 'truth')
    return;
end
for p = 1:numel(names)
    if isfield(data.truth, names{p})
        theta(p) = data.truth.(names{p});
    end
end
end

function names = full_parameter_names()
names = {'Gab','Gbc','Gca','Gammaab','Gammabc','Gammaca','Cab','Cbc','Cca'};
end

function med = median_theta(thetaValues)
med = NaN(1, size(thetaValues, 2));
for idx = 1:size(thetaValues, 2)
    values = thetaValues(:, idx);
    values = values(isfinite(values));
    if ~isempty(values)
        med(idx) = median(values);
    end
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
