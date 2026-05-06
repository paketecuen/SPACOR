function campaign = run_three_wire_blackbox_degradation_campaign(mode)
%RUN_THREE_WIRE_BLACKBOX_DEGRADATION_CAMPAIGN Robustness of 3-wire reports.
%
% Phase 8.2 campaign.  It uses the public three-wire characterization API and
% checks whether the black-box selector remains physically coherent under
% controlled measurement degradations.

if nargin < 1 || isempty(mode)
    mode = 'quick';
end

config = campaign_config();
if strcmpi(mode, 'quick')
    cases = quick_cases();
    degradations = quick_degradations();
    windowCyclesList = 0.5;
else
    cases = full_cases();
    degradations = full_degradations();
    windowCyclesList = [0.25, 0.5, 1.0];
end

outDir = config.resultsDir;
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
figureDir = fullfile(outDir, 'figures', ['three_wire_blackbox_degradation_' char(mode)]);
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

summaryRows = {};
candidateTables = {};
rankingTables = {};
caseIndex = 0;

for c = 1:numel(cases)
    spec = cases(c);
    for d = 1:numel(degradations)
        deg = degradations(d);
        for w = 1:numel(windowCyclesList)
            caseIndex = caseIndex + 1;
            fprintf('3W degradation case %d | %s | %s | %.3g cycles\n', ...
                caseIndex, spec.name, deg.name, windowCyclesList(w));

            data = build_case_data(spec, config, deg, caseIndex);
            opts = struct('candidateMode', 'mixed_blackbox', ...
                'windowCycles', windowCyclesList(w), ...
                'hopCycles', config.hopCycles);
            report = nv3_characterize_three_wire_load(data, config, opts);
            pass = expected_behavior_pass(spec, report);
            bestPower = best_usable_power(report.ranking);

            summaryRows(end + 1, :) = {caseIndex, spec.name, ...
                spec.expectedBehavior, spec.expectedModel, deg.name, ...
                windowCyclesList(w), report.recommendation.status, ...
                report.recommendation.model, report.recommendation.topology, ...
                report.recommendation.terms, report.recommendation.family, ...
                report.recommendation.selectedCoverage, ...
                report.recommendation.adequateCoverage, ...
                report.recommendation.dominantFraction, ...
                report.recommendation.medianSelectedN, ...
                report.recommendation.medianCondition, ...
                report.recommendation.medianSigmaMin, ...
                report.recommendation.medianGeometricVolume, ...
                report.recommendation.medianEquationResidual, ...
                report.recommendation.medianPowerResidual, ...
                report.recommendation.medianEnergyResidual, ...
                report.input.voltageKvlResidual, ...
                report.input.currentKclResidual, bestPower, pass}; %#ok<AGROW>

            candidates = report.candidates;
            candidates.CaseIndex = repmat(caseIndex, height(candidates), 1);
            candidates.TrueCase = repmat({spec.name}, height(candidates), 1);
            candidates.ExpectedBehavior = repmat({spec.expectedBehavior}, ...
                height(candidates), 1);
            candidates.Degradation = repmat({deg.name}, height(candidates), 1);
            candidates.WindowCycles = repmat(windowCyclesList(w), ...
                height(candidates), 1);
            candidateTables{end + 1} = candidates; %#ok<AGROW>

            ranking = report.ranking;
            ranking.CaseIndex = repmat(caseIndex, height(ranking), 1);
            ranking.TrueCase = repmat({spec.name}, height(ranking), 1);
            ranking.ExpectedBehavior = repmat({spec.expectedBehavior}, ...
                height(ranking), 1);
            ranking.Degradation = repmat({deg.name}, height(ranking), 1);
            ranking.WindowCycles = repmat(windowCyclesList(w), ...
                height(ranking), 1);
            rankingTables{end + 1} = ranking; %#ok<AGROW>

            if should_plot(spec, deg, windowCyclesList(w), mode)
                outFile = fullfile(figureDir, sprintf('%03d_%s_%s_w%.2g.png', ...
                    caseIndex, safe_name(spec.name), safe_name(deg.name), ...
                    windowCyclesList(w)));
                nv3_plot_characterization_report(report, data, outFile, ...
                    struct('truthMargin', 0.15, 'maxRanking', 8));
            end
        end
    end
end

summary = cell2table(summaryRows, 'VariableNames', { ...
    'CaseIndex', 'TrueCase', 'ExpectedBehavior', 'ExpectedModel', ...
    'Degradation', 'WindowCycles', 'Status', 'SelectedModel', ...
    'SelectedTopology', 'SelectedTerms', 'SelectedFamily', ...
    'SelectedCoverage', 'AdequateCoverage', 'DominantFraction', ...
    'MedianSelectedN', 'MedianCondition', 'MedianSigmaMin', ...
    'MedianGeometricVolume', 'MedianEquationResidual', ...
    'MedianPowerResidual', 'MedianEnergyResidual', ...
    'VoltageKvlResidual', 'CurrentKclResidual', 'BestUsablePowerResidual', ...
    'Pass'});
candidates = vertcat(candidateTables{:});
ranking = vertcat(rankingTables{:});

prefix = 'three_wire_blackbox_degradation';
if strcmpi(mode, 'quick')
    prefix = 'three_wire_blackbox_degradation_quick';
end
writetable(summary, fullfile(outDir, [prefix '_summary.csv']));
writetable(candidates, fullfile(outDir, [prefix '_candidates.csv']));
writetable(ranking, fullfile(outDir, [prefix '_ranking.csv']));
save(fullfile(outDir, [prefix '.mat']), 'summary', 'candidates', ...
    'ranking', 'config');
plot_degradation_overview(summary, outDir, prefix);

campaign = struct('summary', summary, 'candidates', candidates, ...
    'ranking', ranking, 'config', config, 'outputPrefix', prefix, ...
    'figureDir', figureDir);
disp(summary);
end

function config = campaign_config()
config = nv3_default_config();
config.fs = 20000;
config.f0 = 50;
config.duration = 0.30;
config.windowCycles = 0.5;
config.hopCycles = 0.1;
config.conditionLimit = 1e8;
config.rankTolerance = 1e-8;
config.sigmaFloor = 1e-8;
config.modelResidualFloor = 2.5e-3;
config.energyResidualFloor = 2.5e-3;
config.energyResidualFactor = 5;
config.selectionRelativeTolerance = 0.10;
config.adaptiveHarmonics = true;
config.harmonicSelectionMode = 'residual';
config.harmonicCandidates = [1 3 5 7 9 11 13];
config.activeThreshold = 2e-2;
end

function cases = quick_cases()
cases = [
    make_case('delta GLC full', 'delta_parallel_glc', 'rich_voltage', ...
    'exact_delta', 'delta_glc_full_unbalanced')
    make_case('wye series RLC', 'wye_series_rlc', 'rich_current', ...
    'exact_wye', 'wye_rlgamma_full_unbalanced')
    make_case('delta GLC fundamental', 'delta_parallel_glc', ...
    'fundamental_voltage', 'lower_order_equivalent', 'not_delta_glc')
    make_case('mixed delta+wye', 'mixed_delta_wye', 'rich_current', ...
    'outside_library', 'reject')
    ];
end

function cases = full_cases()
cases = quick_cases();
cases(end + 1) = make_case('delta GL pair', 'delta_parallel_gl', ...
    'rich_voltage', 'exact_delta', 'delta_gl_pair_ab_bc');
cases(end + 1) = make_case('wye series RL', 'wye_series_rl', ...
    'rich_current', 'exact_wye', 'wye_rl_full_unbalanced');
end

function spec = make_case(name, topology, excitation, expectedBehavior, expectedModel)
spec = struct('name', name, 'topology', topology, 'excitation', excitation, ...
    'expectedBehavior', expectedBehavior, 'expectedModel', expectedModel);
end

function degradations = quick_degradations()
degradations = struct('name', {}, 'noiseSnrDb', {}, 'operation', {});
degradations(end + 1) = make_deg('clean', Inf, struct());
degradations(end + 1) = make_deg('noise_60db', 60, struct());
degradations(end + 1) = make_deg('noise_40db', 40, struct());
degradations(end + 1) = make_deg('gain_mismatch_1pct', Inf, ...
    struct('voltageGain', [1.010, 0.990, 1.004], ...
    'currentGain', [0.992, 1.008, 1.000]));
degradations(end + 1) = make_deg('offset_1pct', Inf, ...
    struct('offsetFraction', 0.01));
degradations(end + 1) = make_deg('quantized_12bit', Inf, ...
    struct('quantizationBits', 12));
degradations(end + 1) = make_deg('moving_average_5', Inf, ...
    struct('movingAverageSamples', 5));
end

function degradations = full_degradations()
degradations = quick_degradations();
degradations(end + 1) = make_deg('noise_35db', 35, struct());
degradations(end + 1) = make_deg('offset_2pct', Inf, ...
    struct('offsetFraction', 0.02));
degradations(end + 1) = make_deg('slow_drift_1pct', Inf, ...
    struct('driftFraction', 0.01));
degradations(end + 1) = make_deg('quantized_10bit', Inf, ...
    struct('quantizationBits', 10));
degradations(end + 1) = make_deg('moving_average_11', Inf, ...
    struct('movingAverageSamples', 11));
end

function deg = make_deg(name, noiseSnrDb, operation)
deg = struct('name', name, 'noiseSnrDb', noiseSnrDb, ...
    'operation', operation);
end

function data = build_case_data(spec, config, degradation, caseIndex)
switch spec.topology
    case {'delta_parallel_gl', 'delta_parallel_glc'}
        [vabComponents, vbcComponents] = voltage_components(spec.excitation);
        opts = struct('model', spec.topology, ...
            'truth', delta_truth(spec.expectedModel), ...
            'vabComponents', vabComponents, ...
            'vbcComponents', vbcComponents, ...
            'noiseSnrDb', degradation.noiseSnrDb, ...
            'randomSeed', 5100 + caseIndex);
        data = nv3_synthetic_delta_parallel_case(config, opts);
    case {'wye_series_rl', 'wye_series_rlc'}
        [iaComponents, ibComponents] = current_components(spec.excitation);
        opts = struct('model', spec.topology, ...
            'truth', wye_truth(), ...
            'iaComponents', iaComponents, ...
            'ibComponents', ibComponents, ...
            'noiseSnrDb', degradation.noiseSnrDb, ...
            'randomSeed', 5100 + caseIndex);
        data = nv3_synthetic_wye_series_case(config, opts);
    case 'mixed_delta_wye'
        [iaComponents, ibComponents] = current_components(spec.excitation);
        opts = struct('iaComponents', iaComponents, ...
            'ibComponents', ibComponents, ...
            'noiseSnrDb', degradation.noiseSnrDb, ...
            'randomSeed', 5100 + caseIndex);
        data = nv3_synthetic_mixed_delta_wye_case(config, opts);
    otherwise
        error('Unknown topology "%s".', spec.topology);
end
data.label = [spec.name ' | ' degradation.name];
data.trueCaseName = spec.name;
data.expectedBehavior = spec.expectedBehavior;
data.degradationName = degradation.name;
data = apply_degradation(data, degradation.operation);
end

function truth = delta_truth(expectedModel)
g0 = 1 / 20;
gamma0 = 1 / 0.080;
c0 = 80e-6;
if contains(expectedModel, 'pair_ab_bc')
    g = g0 * [1, 1, 1.45];
    gamma = gamma0 * [1, 1, 0.72];
    c = c0 * [1, 1, 1.35];
else
    g = g0 * [1, 1.55, 0.62];
    gamma = gamma0 * [1, 0.64, 1.75];
    c = c0 * [1, 1.42, 0.58];
end
truth = struct('Gab', g(1), 'Gbc', g(2), 'Gca', g(3), ...
    'Gammaab', gamma(1), 'Gammabc', gamma(2), 'Gammaca', gamma(3), ...
    'Cab', c(1), 'Cbc', c(2), 'Cca', c(3));
truth.Rab = 1 / truth.Gab;
truth.Rbc = 1 / truth.Gbc;
truth.Rca = 1 / truth.Gca;
truth.Lab = 1 / truth.Gammaab;
truth.Lbc = 1 / truth.Gammabc;
truth.Lca = 1 / truth.Gammaca;
end

function truth = wye_truth()
truth = struct('Ra', 0.50, 'Rb', 1.20, 'Rc', 0.80, ...
    'La', 5.0e-3, 'Lb', 8.0e-3, 'Lc', 3.0e-3, ...
    'Gammaa', 1 / 450e-6, 'Gammab', 1 / 320e-6, ...
    'Gammac', 1 / 680e-6);
truth.Ca = 1 / truth.Gammaa;
truth.Cb = 1 / truth.Gammab;
truth.Cc = 1 / truth.Gammac;
end

function [vabComponents, vbcComponents] = voltage_components(excitation)
switch excitation
    case 'fundamental_voltage'
        vabComponents = [100, 1, 0.22];
        vbcComponents = [83, 1, -1.86];
    case 'rich_voltage'
        vabComponents = [100, 1, 0.00; 8, 5, 0.40; 4, 7, -1.00];
        vbcComponents = [92, 1, -2.18; 6, 5, -0.75; 3, 7, 1.15];
    otherwise
        error('Unknown voltage excitation "%s".', excitation);
end
end

function [iaComponents, ibComponents] = current_components(excitation)
switch excitation
    case 'rich_current'
        iaComponents = [10.0, 1, 0.10; 0.90, 5, -0.30; 0.45, 7, 1.20];
        ibComponents = [8.0, 1, -2.05; 0.70, 5, 0.85; 0.25, 7, -1.40];
    otherwise
        error('Unknown current excitation "%s".', excitation);
end
end

function data = apply_degradation(data, operation)
if isempty(fieldnames(operation))
    return;
end
if any(isfield(operation, {'voltageGain', 'currentGain'}))
    data = nv3_apply_measurement_degradation(data, operation);
end
if isfield(operation, 'offsetFraction')
    data = add_offsets(data, operation.offsetFraction);
end
if isfield(operation, 'driftFraction')
    data = add_drift(data, operation.driftFraction);
end
if isfield(operation, 'quantizationBits')
    data = quantize_channels(data, operation.quantizationBits);
end
if isfield(operation, 'movingAverageSamples')
    data = moving_average_channels(data, operation.movingAverageSamples);
end
data.power.measured = three_wire_power(data.vab, data.vbc, data.ia, data.ib);
end

function data = add_offsets(data, fraction)
vOffsets = fraction * rms_many(data.vab, data.vbc, data.vca) * [1, -0.7, 0.4];
iOffsets = fraction * rms_many(data.ia, data.ib, data.ic) * [-0.6, 1, -0.3];
data.vab = data.vab + vOffsets(1);
data.vbc = data.vbc + vOffsets(2);
data.vca = data.vca + vOffsets(3);
data.ia = data.ia + iOffsets(1);
data.ib = data.ib + iOffsets(2);
data.ic = data.ic + iOffsets(3);
end

function data = add_drift(data, fraction)
x = linspace(-0.5, 0.5, numel(data.t)).';
vScale = fraction * rms_many(data.vab, data.vbc, data.vca);
iScale = fraction * rms_many(data.ia, data.ib, data.ic);
data.vab = data.vab + vScale * x;
data.vbc = data.vbc - 0.8 * vScale * x;
data.vca = data.vca + 0.35 * vScale * x;
data.ia = data.ia - 0.6 * iScale * x;
data.ib = data.ib + iScale * x;
data.ic = data.ic - 0.4 * iScale * x;
end

function data = quantize_channels(data, bits)
fields = {'vab', 'vbc', 'vca', 'ia', 'ib', 'ic'};
for idx = 1:numel(fields)
    name = fields{idx};
    x = data.(name);
    fullScale = max(abs(x));
    step = 2 * fullScale / (2 ^ bits - 1);
    data.(name) = step * round(x / step);
end
end

function data = moving_average_channels(data, width)
fields = {'vab', 'vbc', 'vca', 'ia', 'ib', 'ic'};
for idx = 1:numel(fields)
    name = fields{idx};
    data.(name) = movmean(data.(name), width, 'Endpoints', 'shrink');
end
end

function value = rms_many(varargin)
values = [];
for idx = 1:numel(varargin)
    values = [values; varargin{idx}(:)]; %#ok<AGROW>
end
value = sqrt(mean(values .^ 2, 'omitnan'));
end

function p = three_wire_power(vab, vbc, ia, ib)
p = vab(:) .* ia(:) + vbc(:) .* (ia(:) + ib(:));
end

function pass = expected_behavior_pass(spec, report)
coverageOk = report.recommendation.selectedCoverage >= 0.80;
powerGate = max([0.015, 2.5 * report.recommendation.residualGate, ...
    report.recommendation.energyGate]);
if strcmp(report.recommendation.status, 'measurement_inconsistent')
    pass = max(report.recommendation.voltageKvlResidual, ...
        report.recommendation.currentKclResidual) > 5e-3 && ...
        report.recommendation.bestUsablePowerResidual <= 5e-2;
    return;
end
switch spec.expectedBehavior
    case 'exact_delta'
        pass = strcmp(report.recommendation.status, 'accepted') && ...
            coverageOk && strcmp(report.recommendation.model, ...
            spec.expectedModel) && ...
            report.recommendation.medianPowerResidual <= powerGate;
    case 'exact_wye'
        pass = strcmp(report.recommendation.status, 'accepted') && ...
            coverageOk && strcmp(report.recommendation.model, ...
            spec.expectedModel) && ...
            report.recommendation.medianPowerResidual <= powerGate;
    case 'lower_order_equivalent'
        statusOk = any(strcmp(report.recommendation.status, ...
            {'accepted', 'ambiguous_equivalent_family'}));
        pass = statusOk && ...
            coverageOk && ~strcmp(report.recommendation.terms, 'glc') && ...
            report.recommendation.medianPowerResidual <= powerGate;
    case 'outside_library'
        pass = strcmp(report.recommendation.status, 'outside_library') && ...
            report.recommendation.selectedCoverage == 0 && ...
            best_usable_power(report.ranking) > 5e-2;
    otherwise
        pass = false;
end
end

function value = best_usable_power(ranking)
if isempty(ranking) || height(ranking) == 0
    value = Inf;
    return;
end
values = ranking.MedianPowerResidualUsable;
values = values(isfinite(values));
if isempty(values)
    value = Inf;
else
    value = min(values);
end
end

function yes = should_plot(spec, degradation, windowCycles, mode)
if abs(windowCycles - 0.5) > 1e-12
    yes = false;
    return;
end
if strcmpi(mode, 'quick')
    yes = any(strcmp(degradation.name, {'noise_40db', 'offset_1pct', ...
        'moving_average_5'})) || strcmp(spec.expectedBehavior, 'outside_library');
else
    yes = any(strcmp(degradation.name, {'noise_35db', 'offset_2pct', ...
        'moving_average_11'}));
end
end

function name = safe_name(name)
name = regexprep(char(name), '[^A-Za-z0-9]+', '_');
name = regexprep(name, '^_|_$', '');
end

function plot_degradation_overview(summary, outDir, prefix)
cases = unique(summary.TrueCase, 'stable');
degradations = unique(summary.Degradation, 'stable');
coverage = NaN(numel(cases), numel(degradations));
power = NaN(numel(cases), numel(degradations));
pass = NaN(numel(cases), numel(degradations));
for c = 1:numel(cases)
    for d = 1:numel(degradations)
        mask = strcmp(summary.TrueCase, cases{c}) & ...
            strcmp(summary.Degradation, degradations{d}) & ...
            summary.WindowCycles == 0.5;
        if any(mask)
            coverage(c, d) = median(summary.SelectedCoverage(mask), 'omitnan');
            power(c, d) = median(summary.MedianPowerResidual(mask), 'omitnan');
            pass(c, d) = mean(summary.Pass(mask));
        end
    end
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1600 760]);
tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
imagesc(coverage, [0 1]);
colorbar;
title('Selected coverage');
format_axes(cases, degradations);
nexttile;
imagesc(log10(max(power, realmin)));
colorbar;
title('log10 median power residual');
format_axes(cases, degradations);
nexttile;
imagesc(pass, [0 1]);
colorbar;
title('Pass');
format_axes(cases, degradations);
exportgraphics(fig, fullfile(outDir, [prefix '_overview.png']), 'Resolution', 160);
close(fig);
end

function format_axes(cases, degradations)
set(gca, 'XTick', 1:numel(degradations), 'XTickLabel', degradations, ...
    'YTick', 1:numel(cases), 'YTickLabel', cases, ...
    'TickLabelInterpreter', 'none', 'XTickLabelRotation', 55);
end
