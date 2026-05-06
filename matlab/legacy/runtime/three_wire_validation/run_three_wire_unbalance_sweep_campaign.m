function campaign = run_three_wire_unbalance_sweep_campaign(mode)
%RUN_THREE_WIRE_UNBALANCE_SWEEP_CAMPAIGN Load-unbalance and model-order sweep.
%
% The campaign simulates a known delta parallel RLC load and identifies each
% window with a hierarchy of models:
%
%   balanced        : 3 parameters
%   pair_equal_*    : 6 parameters
%   full_unbalanced : 9 parameters
%
% The selected model is the simplest geometrically valid candidate whose
% residual is compatible with the configured practical/noise floor.

if nargin < 1 || isempty(mode)
    mode = 'full';
end

config = nv3_default_config();
config.fs = 20000;
config.f0 = 50;
config.duration = 0.4;
config.conditionLimit = 1e8;
config.adaptiveHarmonics = true;
config.harmonicCandidates = [1 3 5 7 9 11 13];
config.activeThreshold = 1e-3;
config.modelResidualFloor = 2.5e-3;
config.selectionRelativeTolerance = 0.10;

if strcmpi(mode, 'quick')
    excitationNames = {'rich_curve'};
    unbalanceLevels = [0, 0.25, 0.50, 1.00];
    degradationCases = degradation_cases(true);
    windowCyclesList = [0.5, 1.0];
else
    excitationNames = {'weak_curvature', 'rich_curve'};
    unbalanceLevels = [0, 0.05, 0.10, 0.25, 0.50, 1.00];
    degradationCases = degradation_cases(false);
    windowCyclesList = [0.5, 1.0];
end

outDir = config.resultsDir;
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

summaryRows = {};
selectionTables = {};
parameterTables = {};
caseIndex = 0;

for e = 1:numel(excitationNames)
    excitationName = excitationNames{e};
    [vabComponents, vbcComponents] = excitation_components(excitationName);
    for u = 1:numel(unbalanceLevels)
        unbalanceLevel = unbalanceLevels(u);
        truth = unbalance_truth(unbalanceLevel);
        for d = 1:numel(degradationCases)
            degCase = degradationCases(d);
            for wc = 1:numel(windowCyclesList)
                caseIndex = caseIndex + 1;
                fprintf('Unbalance sweep case %d | %s | u=%.2f | %s | %.2g cycles\n', ...
                    caseIndex, excitationName, unbalanceLevel, degCase.name, windowCyclesList(wc));
                data = make_case_data(config, truth, vabComponents, vbcComponents, degCase, caseIndex);
                result = identify_adaptive_order(data, config, degCase, windowCyclesList(wc));
                summaryRows(end + 1, :) = { ...
                    caseIndex, excitationName, unbalanceLevel, degCase.name, ...
                    windowCyclesList(wc), result.windowCount, result.adequateCoverage, ...
                    result.medianSelectedN, result.fracBalanced, result.fracPair, ...
                    result.fracFull, result.fracNoAdequate, result.medianCondition, ...
                    result.medianPowerResidual, result.maxBranchRelativeError, ...
                    result.pass}; %#ok<AGROW>
                selectionTables{end + 1} = result.selectionTable; %#ok<AGROW>
                parameterTables{end + 1} = parameter_summary_table(caseIndex, ...
                    excitationName, unbalanceLevel, degCase.name, windowCyclesList(wc), result); %#ok<AGROW>

                if should_plot_case(excitationName, unbalanceLevel, degCase.name, windowCyclesList(wc))
                    plotFile = sprintf('three_wire_unbalance_%s_u%03d_%s_w%.2g_adaptive.png', ...
                        excitationName, round(100 * unbalanceLevel), degCase.name, windowCyclesList(wc));
                    plot_adaptive_case(data, result, fullfile(outDir, plotFile));
                end
            end
        end
    end
end

summary = cell2table(summaryRows, 'VariableNames', { ...
    'CaseIndex', 'Excitation', 'UnbalanceLevel', 'Degradation', ...
    'WindowCycles', 'WindowCount', 'AdequateCoverage', 'MedianSelectedN', ...
    'FractionBalanced', 'FractionPair', 'FractionFull', 'FractionNoAdequate', ...
    'MedianCondition', 'MedianPowerResidual', 'MaxBranchRelativeError', 'Pass'});
selection = vertcat(selectionTables{:});
parameters = vertcat(parameterTables{:});

prefix = 'three_wire_unbalance_sweep';
if strcmpi(mode, 'quick')
    prefix = 'three_wire_unbalance_sweep_quick';
end
writetable(summary, fullfile(outDir, [prefix '_summary.csv']));
writetable(selection, fullfile(outDir, [prefix '_selection.csv']));
writetable(parameters, fullfile(outDir, [prefix '_parameters.csv']));
save(fullfile(outDir, [prefix '.mat']), 'summary', 'selection', 'parameters');
plot_unbalance_summary(summary, outDir, prefix);

campaign = struct('summary', summary, 'selection', selection, ...
    'parameters', parameters, 'config', config, 'outputPrefix', prefix);
disp(summary);
end

function cases = degradation_cases(quickMode)
cases = struct('name', {}, 'noiseSnrDb', {}, 'degradation', {}, 'quantizationBits', {});
cases(end + 1) = struct('name', 'clean', 'noiseSnrDb', Inf, ...
    'degradation', struct(), 'quantizationBits', Inf);
cases(end + 1) = struct('name', 'noise_60db', 'noiseSnrDb', 60, ...
    'degradation', struct(), 'quantizationBits', Inf);
cases(end + 1) = struct('name', 'noise_40db', 'noiseSnrDb', 40, ...
    'degradation', struct(), 'quantizationBits', Inf);
if ~quickMode
    cases(end + 1) = struct('name', 'gain_mismatch', 'noiseSnrDb', Inf, ...
        'degradation', struct('voltageGain', [1.002 0.998 1.000], ...
        'currentGain', [0.997 1.003 1.000]), 'quantizationBits', Inf);
    cases(end + 1) = struct('name', 'quantized_12bit', 'noiseSnrDb', Inf, ...
        'degradation', struct(), 'quantizationBits', 12);
end
end

function [vabComponents, vbcComponents] = excitation_components(name)
switch name
    case 'weak_curvature'
        vabComponents = [100, 1, 0.00; 0.25, 5, 0.40];
        vbcComponents = [92, 1, -2.18; 0.20, 5, -0.75];
    case 'rich_curve'
        vabComponents = [100, 1, 0.00; 8, 5, 0.40; 4, 7, -1.00];
        vbcComponents = [92, 1, -2.18; 6, 5, -0.75; 3, 7, 1.15];
    otherwise
        error('Unknown excitation "%s".', name);
end
end

function truth = unbalance_truth(level)
g0 = 1 / 20;
gamma0 = 1 / 0.080;
c0 = 80e-6;
gMult = [1, 1 + 0.80 * level, max(0.20, 1 - 0.60 * level)];
gammaMult = [1, max(0.25, 1 - 0.55 * level), 1 + 0.90 * level];
cMult = [1, 1 + 0.70 * level, max(0.25, 1 - 0.50 * level)];
truth = struct();
truth.Gab = g0 * gMult(1);
truth.Gbc = g0 * gMult(2);
truth.Gca = g0 * gMult(3);
truth.Gammaab = gamma0 * gammaMult(1);
truth.Gammabc = gamma0 * gammaMult(2);
truth.Gammaca = gamma0 * gammaMult(3);
truth.Cab = c0 * cMult(1);
truth.Cbc = c0 * cMult(2);
truth.Cca = c0 * cMult(3);
truth.Rab = 1 / truth.Gab;
truth.Rbc = 1 / truth.Gbc;
truth.Rca = 1 / truth.Gca;
truth.Lab = 1 / truth.Gammaab;
truth.Lbc = 1 / truth.Gammabc;
truth.Lca = 1 / truth.Gammaca;
end

function data = make_case_data(config, truth, vabComponents, vbcComponents, degCase, caseIndex)
opts = struct();
opts.model = 'delta_parallel_glc';
opts.truth = truth;
opts.vabComponents = vabComponents;
opts.vbcComponents = vbcComponents;
opts.noiseSnrDb = degCase.noiseSnrDb;
opts.randomSeed = 200 + caseIndex;
data = nv3_synthetic_delta_parallel_case(config, opts);
if isfield(degCase, 'degradation')
    data = nv3_apply_measurement_degradation(data, degCase.degradation);
end
if isfield(degCase, 'quantizationBits') && isfinite(degCase.quantizationBits)
    data = quantize_measurements(data, degCase.quantizationBits);
end
end

function data = quantize_measurements(data, bits)
fields = {'vab', 'vbc', 'vca', 'ia', 'ib', 'ic'};
for idx = 1:numel(fields)
    name = fields{idx};
    x = data.(name);
    fullScale = max(abs(x));
    step = 2 * fullScale / (2 ^ bits - 1);
    data.(name) = step * round(x / step);
end
data.power.measured = data.vab(:) .* data.ia(:) + ...
    data.vbc(:) .* (data.ia(:) + data.ib(:));
data.quantizationBits = bits;
end

function result = identify_adaptive_order(data, config, degCase, windowCycles)
windowSamples = max(12, round(windowCycles * data.fs / data.f0));
hopSamples = max(1, round(0.1 * data.fs / data.f0));
n = numel(data.t);
starts = 1:hopSamples:(n - windowSamples + 1);
candidates = candidate_models();
truthTheta = truth_vector(data.truth);
namesFull = full_parameter_names();
gate = residual_gate(config, degCase);

selectedFullTheta = NaN(numel(starts), 9);
selectionRows = cell(numel(starts), 14);
for w = 1:numel(starts)
    idx = starts(w):(starts(w) + windowSamples - 1);
    signals = fitted_window_signals(data, idx, config);
    [Xfull, y] = full_delta_system(signals);
    evaluations = evaluate_candidates(candidates, Xfull, y, signals, data, idx, config);
    chosen = choose_candidate(evaluations, gate, config.selectionRelativeTolerance);
    if chosen.selected
        selectedFullTheta(w, :) = chosen.fullTheta.';
    end
    selectionRows(w, :) = {w, mean(data.t(idx)), chosen.name, chosen.family, ...
        chosen.nParameters, chosen.rank, chosen.conditionNumber, ...
        chosen.equationResidual, chosen.powerResidual, chosen.adequate, ...
        gate, chosen.passive, chosen.selected, chosen.reason};
end

selectionTable = cell2table(selectionRows, 'VariableNames', { ...
    'Window', 'CenterTime', 'SelectedModel', 'SelectedFamily', ...
    'SelectedN', 'Rank', 'ConditionNumber', 'EquationResidual', ...
    'PowerResidual', 'Adequate', 'ResidualGate', 'Passive', ...
    'Selected', 'Reason'});

adequate = selectionTable.Adequate > 0;
selectedN = selectionTable.SelectedN;
result = struct();
result.selectionTable = selectionTable;
result.selectedFullTheta = selectedFullTheta;
result.windowCount = numel(starts);
result.adequateCoverage = mean(adequate);
result.medianSelectedN = median_finite(selectedN(adequate));
result.fracBalanced = mean(strcmp(selectionTable.SelectedFamily(adequate), 'balanced'));
result.fracPair = mean(strcmp(selectionTable.SelectedFamily(adequate), 'pair_equal'));
result.fracFull = mean(strcmp(selectionTable.SelectedFamily(adequate), 'full_unbalanced'));
result.fracNoAdequate = mean(~adequate);
result.medianCondition = median_finite(selectionTable.ConditionNumber(adequate));
result.medianPowerResidual = median_finite(selectionTable.PowerResidual(adequate));
result.truthTheta = truthTheta;
result.parameterNames = namesFull;
result.parameterMedian = median_theta(selectedFullTheta(adequate, :));
result.parameterRelativeError = abs(result.parameterMedian(:) - truthTheta(:)) ./ max(abs(truthTheta(:)), eps);
result.maxBranchRelativeError = max(result.parameterRelativeError, [], 'omitnan');
result.pass = result.adequateCoverage >= 0.80 && result.maxBranchRelativeError <= 0.05;
end

function gate = residual_gate(config, degCase)
noiseGate = 0;
if isfinite(degCase.noiseSnrDb)
    noiseGate = 5 * 10 ^ (-degCase.noiseSnrDb / 20);
end
if strcmp(degCase.name, 'quantized_12bit')
    noiseGate = max(noiseGate, 2e-3);
elseif strcmp(degCase.name, 'gain_mismatch')
    noiseGate = max(noiseGate, 5e-3);
end
gate = max(config.modelResidualFloor, noiseGate);
end

function signals = fitted_window_signals(data, idx, config)
t = data.t(idx);
raw = [data.vab(idx), data.vbc(idx), data.vca(idx), ...
    data.ia(idx), data.ib(idx), data.ic(idx)];
selection = nv3_select_active_harmonics(raw, t, data.f0, config.harmonicCandidates, config);
h = selection.selectedHarmonics;
vabFit = nv3_harmonic_fit(data.vab(idx), t, data.f0, h, struct('maxOrder', 1));
vbcFit = nv3_harmonic_fit(data.vbc(idx), t, data.f0, h, struct('maxOrder', 1));
iaFit = nv3_harmonic_fit(data.ia(idx), t, data.f0, h, struct('maxOrder', 0));
ibFit = nv3_harmonic_fit(data.ib(idx), t, data.f0, h, struct('maxOrder', 0));
icFit = nv3_harmonic_fit(data.ic(idx), t, data.f0, h, struct('maxOrder', 0));
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
signals.harmonics = h;
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

function candidates = candidate_models()
candidates = struct('name', {}, 'family', {}, 'S', {}, 'nParameters', {});
candidates(end + 1) = make_candidate('balanced', 'balanced', symmetry_matrix('balanced'));
candidates(end + 1) = make_candidate('pair_ab_bc', 'pair_equal', symmetry_matrix('pair_ab_bc'));
candidates(end + 1) = make_candidate('pair_bc_ca', 'pair_equal', symmetry_matrix('pair_bc_ca'));
candidates(end + 1) = make_candidate('pair_ca_ab', 'pair_equal', symmetry_matrix('pair_ca_ab'));
candidates(end + 1) = make_candidate('full_unbalanced', 'full_unbalanced', eye(9));
end

function candidate = make_candidate(name, family, S)
candidate = struct('name', name, 'family', family, 'S', S, 'nParameters', size(S, 2));
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

function evaluations = evaluate_candidates(candidates, Xfull, y, signals, data, idx, config)
evaluations = candidates;
for k = 1:numel(candidates)
    X = Xfull * candidates(k).S;
    [theta, fit] = normalized_solve(X, y);
    fullTheta = candidates(k).S * theta;
    passive = all(isfinite(theta)) && all(theta > 0);
    rankOk = fit.rank == candidates(k).nParameters;
    conditionOk = isfinite(fit.conditionNumber) && fit.conditionNumber <= config.conditionLimit;
    powerResidual = Inf;
    if rankOk && conditionOk && passive
        predPower = predict_power(signals, fullTheta);
        measuredPower = data.vab(idx) .* data.ia(idx) + data.vbc(idx) .* ...
            (data.ia(idx) + data.ib(idx));
        powerResidual = norm(measuredPower(:) - predPower(:)) / max(norm(measuredPower(:)), eps);
    end
    evaluations(k).theta = theta;
    evaluations(k).fullTheta = fullTheta;
    evaluations(k).rank = fit.rank;
    evaluations(k).conditionNumber = fit.conditionNumber;
    evaluations(k).equationResidual = fit.relativeResidual;
    evaluations(k).powerResidual = powerResidual;
    evaluations(k).passive = passive;
    evaluations(k).geometricValid = rankOk && conditionOk;
    evaluations(k).usable = rankOk && conditionOk && passive && isfinite(powerResidual);
end
end

function chosen = choose_candidate(evaluations, gate, relativeTolerance)
families = {'balanced', 'pair_equal', 'full_unbalanced'};
chosen = empty_choice('no_geometric_candidate');
usableAll = evaluations([evaluations.usable]);
if isempty(usableAll)
    return;
end
allScores = [usableAll.powerResidual] + [usableAll.equationResidual];
bestScore = min(allScores);
for f = 1:numel(families)
    family = families{f};
    group = evaluations(strcmp({evaluations.family}, family));
    group = group([group.usable]);
    if isempty(group)
        continue;
    end
    [~, bestIdx] = min([group.powerResidual] + [group.equationResidual]);
    candidate = group(bestIdx);
    candidateScore = candidate.powerResidual + candidate.equationResidual;
    scoreClose = candidateScore <= (1 + relativeTolerance) * bestScore + 0.10 * gate;
    adequate = candidate.powerResidual <= gate && candidate.equationResidual <= 2 * gate && scoreClose;
    if adequate
        chosen = candidate_to_choice(candidate, true, 'adequate_smallest_model');
        return;
    end
    if ~chosen.selected
        chosen = candidate_to_choice(candidate, false, 'best_not_adequate');
    end
end
end

function choice = empty_choice(reason)
choice = struct('name', 'none', 'family', 'none', 'nParameters', NaN, ...
    'rank', NaN, 'conditionNumber', NaN, 'equationResidual', NaN, ...
    'powerResidual', NaN, 'adequate', false, 'passive', false, ...
    'selected', false, 'reason', reason, 'fullTheta', NaN(9, 1));
end

function choice = candidate_to_choice(candidate, adequate, reason)
choice = struct('name', candidate.name, 'family', candidate.family, ...
    'nParameters', candidate.nParameters, 'rank', candidate.rank, ...
    'conditionNumber', candidate.conditionNumber, ...
    'equationResidual', candidate.equationResidual, ...
    'powerResidual', candidate.powerResidual, 'adequate', adequate, ...
    'passive', candidate.passive, 'selected', true, 'reason', reason, ...
    'fullTheta', candidate.fullTheta);
end

function [theta, fit] = normalized_solve(X, y)
mask = all(isfinite(X), 2) & isfinite(y);
Xv = X(mask, :);
yv = y(mask);
m = size(X, 2);
theta = NaN(m, 1);
fit = struct('rank', 0, 'conditionNumber', Inf, 'relativeResidual', Inf);
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

function p = predict_power(s, theta)
iab = theta(1) * s.vab + theta(4) * s.qab + theta(7) * s.dvab;
ibc = theta(2) * s.vbc + theta(5) * s.qbc + theta(8) * s.dvbc;
ica = theta(3) * s.vca + theta(6) * s.qca + theta(9) * s.dvca;
ia = iab - ica;
ib = ibc - iab;
p = s.vab(:) .* ia(:) + s.vbc(:) .* (ia(:) + ib(:));
end

function theta = truth_vector(truth)
theta = [truth.Gab; truth.Gbc; truth.Gca; ...
    truth.Gammaab; truth.Gammabc; truth.Gammaca; ...
    truth.Cab; truth.Cbc; truth.Cca];
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

function tableOut = parameter_summary_table(caseIndex, excitationName, unbalanceLevel, degradation, windowCycles, result)
rows = cell(numel(result.parameterNames), 9);
for p = 1:numel(result.parameterNames)
    rows(p, :) = {caseIndex, excitationName, unbalanceLevel, degradation, ...
        windowCycles, result.parameterNames{p}, result.truthTheta(p), ...
        result.parameterMedian(p), result.parameterRelativeError(p)};
end
tableOut = cell2table(rows, 'VariableNames', { ...
    'CaseIndex', 'Excitation', 'UnbalanceLevel', 'Degradation', ...
    'WindowCycles', 'Parameter', 'Truth', 'MedianEstimate', 'RelativeError'});
end

function yes = should_plot_case(excitationName, unbalanceLevel, degradation, windowCycles)
yes = strcmp(excitationName, 'rich_curve') && abs(unbalanceLevel - 0.5) < 1e-12 && ...
    any(strcmp(degradation, {'clean', 'noise_40db'})) && abs(windowCycles - 0.5) < 1e-12;
end

function plot_adaptive_case(data, result, outFile)
windows = result.selectionTable;
valid = windows.Adequate > 0;
theta = result.selectedFullTheta;
truth = result.truthTheta;
names = result.parameterNames;
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1600 1200]);
tiledlayout(fig, 5, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile([1 2]);
plot(data.t, data.vab, 'LineWidth', 0.8);
hold on;
plot(data.t, data.vbc, 'LineWidth', 0.8);
plot(data.t, data.vca, 'LineWidth', 0.8);
grid on;
title('Line voltages');
legend({'v_{ab}','v_{bc}','v_{ca}'}, 'Location', 'best');
nexttile;
plot(data.t, data.ia, 'LineWidth', 0.8);
hold on;
plot(data.t, data.ib, 'LineWidth', 0.8);
plot(data.t, data.ic, 'LineWidth', 0.8);
grid on;
title('Line currents');
legend({'i_a','i_b','i_c'}, 'Location', 'best');
nexttile([1 3]);
stairs(windows.CenterTime, windows.SelectedN, 'LineWidth', 1.4);
hold on;
plot(windows.CenterTime(~valid), windows.SelectedN(~valid), 'rx');
grid on;
ylim([0 10]);
ylabel('selected n_\theta');
title('Adaptive parameter count per window');
for p = 1:numel(names)
    nexttile;
    plot(windows.CenterTime(valid), theta(valid, p), '-o', 'MarkerSize', 3, 'LineWidth', 1.0);
    hold on;
    yline(truth(p), '--', 'Truth');
    grid on;
    title(names{p}, 'Interpreter', 'none');
    apply_truth_zoom(theta(valid, p), truth(p), 0.10);
end
exportgraphics(fig, outFile, 'Resolution', 160);
close(fig);
end

function plot_unbalance_summary(summary, outDir, prefix)
excitations = unique(summary.Excitation, 'stable');
degradations = unique(summary.Degradation, 'stable');
for e = 1:numel(excitations)
    for d = 1:numel(degradations)
        mask = strcmp(summary.Excitation, excitations{e}) & ...
            strcmp(summary.Degradation, degradations{d}) & summary.WindowCycles == 0.5;
        if ~any(mask)
            continue;
        end
        sub = summary(mask, :);
        [levels, order] = sort(sub.UnbalanceLevel);
        sub = sub(order, :);
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 800]);
        tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
        nexttile;
        plot(levels, sub.AdequateCoverage, '-o', 'LineWidth', 1.3);
        grid on;
        ylim([0 1.05]);
        xlabel('Load unbalance level');
        ylabel('adequate coverage');
        title('Coverage');
        nexttile;
        plot(levels, sub.MedianSelectedN, '-o', 'LineWidth', 1.3);
        grid on;
        ylim([0 10]);
        xlabel('Load unbalance level');
        ylabel('median selected n_\theta');
        title('Selected model order');
        nexttile;
        plot(levels, sub.FractionBalanced, '-o', levels, sub.FractionPair, '-o', ...
            levels, sub.FractionFull, '-o', 'LineWidth', 1.3);
        grid on;
        ylim([0 1.05]);
        xlabel('Load unbalance level');
        ylabel('fraction');
        legend({'balanced','pair','full'}, 'Location', 'best');
        title('Selected families');
        nexttile;
        semilogy(levels, sub.MaxBranchRelativeError, '-o', 'LineWidth', 1.3);
        grid on;
        xlabel('Load unbalance level');
        ylabel('max branch relative error');
        title('Selected branch-parameter error');
        sgtitle(sprintf('%s | %s | 0.5 cycles', excitations{e}, degradations{d}), ...
            'Interpreter', 'none');
        exportgraphics(fig, fullfile(outDir, sprintf('%s_%s_%s_w0.5.png', ...
            prefix, excitations{e}, degradations{d})), 'Resolution', 160);
        close(fig);
    end
end
end

function apply_truth_zoom(values, truthValue, margin)
values = values(isfinite(values));
if isempty(values) || ~isfinite(truthValue)
    return;
end
p10 = percentile(values, 10);
p90 = percentile(values, 90);
dev = max([abs(p10 - truthValue), abs(p90 - truthValue), ...
    margin * max(abs(truthValue), eps), eps]);
ylim([truthValue - 1.15 * dev, truthValue + 1.15 * dev]);
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
