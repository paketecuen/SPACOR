function campaign = run_three_wire_blackbox_equivalent_campaign(mode)
%RUN_THREE_WIRE_BLACKBOX_EQUIVALENT_CAMPAIGN Black-box 3-wire equivalents.
%
% The simulated load is known, but the identifier is only allowed to see
% terminal voltages/currents.  It chooses among delta parallel equivalents:
%
%   G, Gamma, C, G+Gamma, G+C, Gamma+C, G+Gamma+C
%
% combined with balanced, pair-equal and full-unbalanced branch symmetries,
% and among wye series R/L equivalents with the same symmetry hierarchy.

if nargin < 1 || isempty(mode)
    mode = 'quick';
end

config = nv3_default_config();
config.fs = 20000;
config.f0 = 50;
config.duration = 0.4;
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

if strcmpi(mode, 'quick')
    cases = quick_cases();
    degradationNames = {'clean', 'noise_60db'};
else
    cases = full_cases();
    degradationNames = {'clean', 'noise_60db', 'noise_40db'};
end

outDir = config.resultsDir;
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

summaryRows = {};
selectionTables = {};
candidateTables = {};
rankingTables = {};
caseIndex = 0;

for c = 1:numel(cases)
    spec = cases(c);
    for d = 1:numel(degradationNames)
        caseIndex = caseIndex + 1;
        degradationName = degradationNames{d};
        fprintf('Black-box case %d | %s | %s | %s\n', ...
            caseIndex, spec.name, spec.excitation, degradationName);

        data = build_case_data(spec, config, degradationName, caseIndex);
        result = nv3_identify_adaptive_geometric(data, config, ...
            struct('candidateMode', 'mixed_blackbox', ...
            'windowCycles', config.windowCycles));

        dominantModel = dominant_selected(result.windows, 'SelectedModel');
        dominantTopology = dominant_selected(result.windows, 'SelectedTopology');
        dominantTerms = dominant_selected(result.windows, 'SelectedTerms');
        dominantFamily = dominant_selected(result.windows, 'SelectedFamily');
        dominantFraction = dominant_fraction(result.windows, dominantModel);
        pass = blackbox_pass(spec, result, dominantModel, dominantTerms);

        summaryRows(end + 1, :) = { ...
            caseIndex, spec.name, spec.trueTopology, spec.excitation, ...
            degradationName, spec.expectedBehavior, spec.expectedModel, ...
            height(result.windows), result.summary.SelectedCoverage, ...
            result.summary.AdequateCoverage, result.summary.MedianSelectedN, ...
            dominantModel, dominantTopology, dominantTerms, dominantFamily, ...
            dominantFraction, result.summary.MedianCondition, result.summary.MedianSigmaMin, ...
            result.summary.MedianGeometricVolume, ...
            result.summary.MedianEquationResidual, ...
            result.summary.MedianPowerResidual, ...
            result.summary.MedianEnergyResidual, pass}; %#ok<AGROW>

        selection = result.windows;
        selection.CaseIndex = repmat(caseIndex, height(selection), 1);
        selection.TrueCase = repmat({spec.name}, height(selection), 1);
        selection.TrueTopology = repmat({spec.trueTopology}, height(selection), 1);
        selection.Excitation = repmat({spec.excitation}, height(selection), 1);
        selection.Degradation = repmat({degradationName}, height(selection), 1);
        selectionTables{end + 1} = selection; %#ok<AGROW>

        candidates = result.candidates;
        candidates.CaseIndex = repmat(caseIndex, height(candidates), 1);
        candidates.TrueCase = repmat({spec.name}, height(candidates), 1);
        candidates.TrueTopology = repmat({spec.trueTopology}, height(candidates), 1);
        candidates.Excitation = repmat({spec.excitation}, height(candidates), 1);
        candidates.Degradation = repmat({degradationName}, height(candidates), 1);
        candidateTables{end + 1} = candidates; %#ok<AGROW>

        rankingTables{end + 1} = candidate_ranking(caseIndex, spec, ...
            degradationName, result.candidates); %#ok<AGROW>
    end
end

summary = cell2table(summaryRows, 'VariableNames', { ...
    'CaseIndex', 'TrueCase', 'TrueTopology', 'Excitation', ...
    'Degradation', 'ExpectedBehavior', 'ExpectedModel', 'WindowCount', ...
    'SelectedCoverage', 'AdequateCoverage', 'MedianSelectedN', ...
    'DominantModel', 'DominantTopology', 'DominantTerms', 'DominantFamily', ...
    'DominantModelFraction', 'MedianCondition', 'MedianSigmaMin', ...
    'MedianGeometricVolume', 'MedianEquationResidual', ...
    'MedianPowerResidual', 'MedianEnergyResidual', 'Pass'});
selection = vertcat(selectionTables{:});
candidates = vertcat(candidateTables{:});
ranking = vertcat(rankingTables{:});

prefix = 'three_wire_blackbox_equivalent';
if strcmpi(mode, 'quick')
    prefix = 'three_wire_blackbox_equivalent_quick';
end
writetable(summary, fullfile(outDir, [prefix '_summary.csv']));
writetable(selection, fullfile(outDir, [prefix '_selection.csv']));
writetable(candidates, fullfile(outDir, [prefix '_candidates.csv']));
writetable(ranking, fullfile(outDir, [prefix '_ranking.csv']));
save(fullfile(outDir, [prefix '.mat']), 'summary', 'selection', ...
    'candidates', 'ranking', 'config');
plot_blackbox_summary(summary, outDir, prefix);

campaign = struct('summary', summary, 'selection', selection, ...
    'candidates', candidates, 'ranking', ranking, 'config', config, ...
    'outputPrefix', prefix);
disp(summary);
end

function cases = quick_cases()
cases = [
    make_case('delta G balanced', 'delta_parallel_g', 'rich_curve', ...
    'terminal_equivalent', 'delta_g_balanced|wye_r_balanced')
    make_case('delta GL pair', 'delta_parallel_gl', 'rich_curve', ...
    'exact_delta', 'delta_gl_pair_ab_bc')
    make_case('delta GLC full', 'delta_parallel_glc', 'rich_curve', ...
    'exact_delta', 'delta_glc_full_unbalanced')
    make_case('delta GLC fundamental', 'delta_parallel_glc', ...
    'fundamental_unbalanced', 'lower_order_equivalent', 'not_delta_glc')
    make_case('wye series RL', 'wye_series_rl', 'rich_curve_current', ...
    'exact_wye', 'wye_rl_full_unbalanced')
    make_case('wye series RLC', 'wye_series_rlc', 'rich_curve_current', ...
    'exact_wye', 'wye_rlgamma_full_unbalanced')
    make_case('mixed delta+wye', 'mixed_delta_wye', 'rich_curve_current', ...
    'outside_library', 'no_pure_delta_or_wye_equivalent')
    ];
end

function cases = full_cases()
cases = quick_cases();
cases(end + 1) = make_case('delta GC full', 'delta_parallel_gc', ...
    'rich_curve', 'exact_delta', 'delta_gc_full_unbalanced');
cases(end + 1) = make_case('delta GammaC full', 'delta_parallel_gammac', ...
    'rich_curve', 'exact_delta', 'delta_gammac_full_unbalanced');
end

function spec = make_case(name, trueTopology, excitation, expectedBehavior, expectedModel)
spec = struct('name', name, 'trueTopology', trueTopology, ...
    'excitation', excitation, 'expectedBehavior', expectedBehavior, ...
    'expectedModel', expectedModel);
end

function data = build_case_data(spec, config, degradationName, caseIndex)
switch degradationName
    case 'clean'
        noiseSnrDb = Inf;
    case 'noise_60db'
        noiseSnrDb = 60;
    case 'noise_40db'
        noiseSnrDb = 40;
    otherwise
        error('Unknown degradation "%s".', degradationName);
end

switch spec.trueTopology
    case {'delta_parallel_g', 'delta_parallel_gl', 'delta_parallel_gc', ...
            'delta_parallel_gammac', 'delta_parallel_glc'}
        [vabComponents, vbcComponents] = voltage_components(spec.excitation);
        opts = struct('model', spec.trueTopology, ...
            'truth', delta_truth(spec.expectedModel), ...
            'vabComponents', vabComponents, ...
            'vbcComponents', vbcComponents, ...
            'noiseSnrDb', noiseSnrDb, ...
            'randomSeed', 1200 + caseIndex);
        data = nv3_synthetic_delta_parallel_case(config, opts);
    case {'wye_series_rl', 'wye_series_rlc'}
        [iaComponents, ibComponents] = current_components(spec.excitation);
        opts = struct('model', spec.trueTopology, ...
            'truth', wye_truth(), ...
            'iaComponents', iaComponents, ...
            'ibComponents', ibComponents, ...
            'noiseSnrDb', noiseSnrDb, ...
            'randomSeed', 1200 + caseIndex);
        data = nv3_synthetic_wye_series_case(config, opts);
    case 'mixed_delta_wye'
        [iaComponents, ibComponents] = current_components(spec.excitation);
        opts = struct('iaComponents', iaComponents, ...
            'ibComponents', ibComponents, ...
            'noiseSnrDb', noiseSnrDb, ...
            'randomSeed', 1200 + caseIndex);
        data = nv3_synthetic_mixed_delta_wye_case(config, opts);
    otherwise
        error('Unknown true topology "%s".', spec.trueTopology);
end
data.trueCaseName = spec.name;
data.trueTopology = spec.trueTopology;
data.blackboxExcitation = spec.excitation;
data.degradationName = degradationName;
end

function truth = delta_truth(expectedModel)
g0 = 1 / 20;
gamma0 = 1 / 0.080;
c0 = 80e-6;
if contains(expectedModel, 'full_unbalanced') || strcmp(expectedModel, 'not_delta_glc')
    g = g0 * [1, 1.55, 0.62];
    gamma = gamma0 * [1, 0.64, 1.75];
    c = c0 * [1, 1.42, 0.58];
elseif contains(expectedModel, 'pair_ab_bc')
    g = g0 * [1, 1, 1.45];
    gamma = gamma0 * [1, 1, 0.72];
    c = c0 * [1, 1, 1.35];
elseif contains(expectedModel, 'balanced')
    g = g0 * [1, 1, 1];
    gamma = gamma0 * [1, 1, 1];
    c = c0 * [1, 1, 1];
else
    error('Unknown expected model "%s".', expectedModel);
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
    case 'fundamental_unbalanced'
        vabComponents = [100, 1, 0.22];
        vbcComponents = [83, 1, -1.86];
    case 'rich_curve'
        vabComponents = [100, 1, 0.00; 8, 5, 0.40; 4, 7, -1.00];
        vbcComponents = [92, 1, -2.18; 6, 5, -0.75; 3, 7, 1.15];
    otherwise
        error('Unknown voltage excitation "%s".', excitation);
end
end

function [iaComponents, ibComponents] = current_components(excitation)
switch excitation
    case 'rich_curve_current'
        iaComponents = [10.0, 1, 0.10; 0.90, 5, -0.30; 0.45, 7, 1.20];
        ibComponents = [8.0, 1, -2.05; 0.70, 5, 0.85; 0.25, 7, -1.40];
    otherwise
        error('Unknown current excitation "%s".', excitation);
end
end

function pass = blackbox_pass(spec, result, dominantModel, dominantTerms)
switch spec.expectedBehavior
    case 'terminal_equivalent'
        allowedModels = strsplit(spec.expectedModel, '|');
        pass = result.summary.SelectedCoverage >= 0.80 && ...
            any(strcmp(dominantModel, allowedModels)) && ...
            result.summary.MedianPowerResidual <= 5e-3;
    case 'exact_delta'
        pass = result.summary.SelectedCoverage >= 0.80 && ...
            strcmp(dominantModel, spec.expectedModel) && ...
            result.summary.MedianPowerResidual <= 5e-3;
    case 'exact_wye'
        pass = result.summary.SelectedCoverage >= 0.80 && ...
            strcmp(dominantModel, spec.expectedModel) && ...
            result.summary.MedianPowerResidual <= 5e-3;
    case 'lower_order_equivalent'
        pass = result.summary.SelectedCoverage >= 0.80 && ...
            ~strcmp(dominantTerms, 'glc') && ...
            result.summary.MedianPowerResidual <= 5e-3;
    case 'outside_library'
        pass = result.summary.SelectedCoverage == 0 && ...
            best_usable_power_residual(result.candidates) > 5e-2;
    case 'topology_mismatch'
        pass = result.summary.SelectedCoverage == 0 || ...
            result.summary.MedianPowerResidual <= 5e-2;
    otherwise
        pass = false;
end
end

function value = best_usable_power_residual(candidates)
usable = candidates.Usable > 0 & isfinite(candidates.PowerResidual);
if ~any(usable)
    value = Inf;
else
    value = min(candidates.PowerResidual(usable));
end
end

function value = dominant_selected(windows, fieldName)
selected = windows.Selected > 0;
if ~any(selected)
    value = 'reject';
    return;
end
values = windows.(fieldName)(selected);
u = unique(values, 'stable');
counts = zeros(numel(u), 1);
for idx = 1:numel(u)
    counts(idx) = sum(strcmp(values, u{idx}));
end
[~, maxIdx] = max(counts);
value = u{maxIdx};
end

function fraction = dominant_fraction(windows, modelName)
if strcmp(modelName, 'reject')
    fraction = mean(windows.Selected == 0);
    return;
end
selected = windows.Selected > 0;
if ~any(selected)
    fraction = 0;
else
    fraction = mean(strcmp(windows.SelectedModel(selected), modelName));
end
end

function ranking = candidate_ranking(caseIndex, spec, degradationName, candidates)
models = unique(candidates.CandidateModel, 'stable');
rows = cell(numel(models), 14);
for m = 1:numel(models)
    model = models{m};
    mask = strcmp(candidates.CandidateModel, model);
    sub = candidates(mask, :);
    rows(m, :) = {caseIndex, spec.name, spec.trueTopology, spec.excitation, ...
        degradationName, model, sub.CandidateTopology{1}, ...
        sub.CandidateTerms{1}, sub.CandidateFamily{1}, sub.NParameters(1), ...
        mean(sub.Usable > 0), ...
        median_finite(sub.EquationResidual(sub.Usable > 0)), ...
        median_finite(sub.PowerResidual(sub.Usable > 0)), ...
        median_finite(sub.EnergyResidual(sub.Usable > 0))};
end
ranking = cell2table(rows, 'VariableNames', { ...
    'CaseIndex', 'TrueCase', 'TrueTopology', 'Excitation', ...
    'Degradation', 'CandidateModel', 'CandidateTopology', 'CandidateTerms', ...
    'CandidateFamily', 'NParameters', 'UsableCoverage', ...
    'MedianEquationResidual', 'MedianPowerResidual', ...
    'MedianEnergyResidual'});
end

function plot_blackbox_summary(summary, outDir, prefix)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 900]);
tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
x = 1:height(summary);
labels = compose('%d %s %s', summary.CaseIndex, ...
    string(summary.TrueCase), string(summary.Degradation));

nexttile;
bar(x, summary.SelectedCoverage);
grid on;
ylim([0 1.05]);
ylabel('selected coverage');
title('Accepted windows');

nexttile;
bar(x, summary.MedianSelectedN);
grid on;
ylim([0 10]);
ylabel('median n_\theta');
title('Equivalent order');

nexttile;
semilogy(x, summary.MedianPowerResidual, 'o-', 'LineWidth', 1.1);
grid on;
ylabel('median power residual');
title('Terminal power consistency');

nexttile;
bar(x, double(summary.Pass));
grid on;
ylim([0 1.05]);
ylabel('pass');
title('Expected behavior check');

for ax = findall(fig, 'Type', 'axes').'
    ax.XTick = x;
    ax.XTickLabel = labels;
    ax.XTickLabelRotation = 65;
    ax.FontSize = 8;
end
exportgraphics(fig, fullfile(outDir, [prefix '_summary.png']), 'Resolution', 160);
close(fig);
end

function value = median_finite(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = median(x);
end
end
