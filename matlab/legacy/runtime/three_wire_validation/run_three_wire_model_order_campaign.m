function campaign = run_three_wire_model_order_campaign(mode)
%RUN_THREE_WIRE_MODEL_ORDER_CAMPAIGN Validate adaptive geometric model order.
%
% The campaign checks whether the geometric selector chooses the right model
% complexity for known delta parallel RLC loads:
%
%   reject          : geometrically poor RLC window
%   balanced        : 3 free parameters
%   pair_equal      : 6 free parameters
%   full_unbalanced : 9 free parameters

if nargin < 1 || isempty(mode)
    mode = 'quick';
end

config = nv3_default_config();
config.fs = 20000;
config.f0 = 50;
config.duration = 0.4;
config.conditionLimit = 1e8;
config.rankTolerance = 1e-8;
config.sigmaFloor = 1e-8;
config.modelResidualFloor = 2.5e-3;
config.energyResidualFloor = 2.5e-3;
config.energyResidualFactor = 5;
config.selectionRelativeTolerance = 0.10;
config.adaptiveHarmonics = true;
config.harmonicCandidates = [1 3 5 7 9 11 13];
% This campaign explicitly tests that noise must not "repair" a planar
% fundamental trajectory.  Use a conservative basis-energy gate so random
% noise harmonics are not promoted into geometric curvature.
config.activeThreshold = 2e-2;
config.harmonicSelectionMode = 'residual';
config.maxBasisCondition = 1e10;
config.windowCycles = 0.5;
config.hopCycles = 0.1;

if strcmpi(mode, 'quick')
    excitationNames = {'fundamental_unbalanced', 'rich_curve'};
    truthNames = {'balanced', 'pair_ab_bc', 'full_unbalanced'};
    degradationCases = {'clean', 'noise_60db'};
    windowCyclesList = 0.5;
else
    excitationNames = {'fundamental_balanced', 'fundamental_unbalanced', ...
        'weak_curvature', 'rich_curve'};
    truthNames = {'balanced', 'pair_ab_bc', 'full_unbalanced'};
    degradationCases = {'clean', 'noise_60db', 'noise_40db'};
    windowCyclesList = [0.5, 1.0];
end

outDir = config.resultsDir;
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

summaryRows = {};
selectionTables = {};
candidateTables = {};
parameterTables = {};
caseIndex = 0;

for e = 1:numel(excitationNames)
    excitationName = excitationNames{e};
    [vabComponents, vbcComponents] = excitation_components(excitationName);
    for tIdx = 1:numel(truthNames)
        truthName = truthNames{tIdx};
        truth = truth_for_family(truthName);
        for d = 1:numel(degradationCases)
            degradationName = degradationCases{d};
            for wc = 1:numel(windowCyclesList)
                caseIndex = caseIndex + 1;
                windowCycles = windowCyclesList(wc);
                fprintf('Model-order case %d | %s | %s | %s | %.2g cycles\n', ...
                    caseIndex, excitationName, truthName, degradationName, windowCycles);

                data = make_case(config, truth, vabComponents, vbcComponents, ...
                    degradationName, caseIndex);
                result = nv3_identify_adaptive_geometric(data, config, ...
                    struct('windowCycles', windowCycles));

                expectedFamily = expected_family(excitationName, truthName);
                dominantFamily = dominant_selected_family(result.windows);
                dominantFraction = dominant_selected_fraction(result.windows, dominantFamily);
                pass = case_passes(result, expectedFamily, dominantFamily, dominantFraction);

                summaryRows(end + 1, :) = { ...
                    caseIndex, excitationName, truthName, degradationName, ...
                    windowCycles, expectedFamily, height(result.windows), ...
                    result.summary.SelectedCoverage, result.summary.AdequateCoverage, ...
                    result.summary.MedianSelectedN, dominantFamily, dominantFraction, ...
                    result.summary.MedianCondition, result.summary.MedianSigmaMin, ...
                    result.summary.MedianGeometricVolume, ...
                    result.summary.MedianEquationResidual, ...
                    result.summary.MedianPowerResidual, ...
                    result.summary.MedianEnergyResidual, ...
                    result.summary.MaxBranchRelativeError, pass}; %#ok<AGROW>

                selection = result.windows;
                selection.CaseIndex = repmat(caseIndex, height(selection), 1);
                selection.Excitation = repmat({excitationName}, height(selection), 1);
                selection.TruthFamily = repmat({truthName}, height(selection), 1);
                selection.Degradation = repmat({degradationName}, height(selection), 1);
                selection.WindowCycles = repmat(windowCycles, height(selection), 1);
                selectionTables{end + 1} = selection; %#ok<AGROW>

                candidates = result.candidates;
                candidates.CaseIndex = repmat(caseIndex, height(candidates), 1);
                candidates.Excitation = repmat({excitationName}, height(candidates), 1);
                candidates.TruthFamily = repmat({truthName}, height(candidates), 1);
                candidates.Degradation = repmat({degradationName}, height(candidates), 1);
                candidates.WindowCycles = repmat(windowCycles, height(candidates), 1);
                candidateTables{end + 1} = candidates; %#ok<AGROW>

                params = result.parameterSummary;
                params.CaseIndex = repmat(caseIndex, height(params), 1);
                params.Excitation = repmat({excitationName}, height(params), 1);
                params.TruthFamily = repmat({truthName}, height(params), 1);
                params.Degradation = repmat({degradationName}, height(params), 1);
                params.WindowCycles = repmat(windowCycles, height(params), 1);
                parameterTables{end + 1} = params; %#ok<AGROW>
            end
        end
    end
end

summary = cell2table(summaryRows, 'VariableNames', { ...
    'CaseIndex', 'Excitation', 'TruthFamily', 'Degradation', ...
    'WindowCycles', 'ExpectedFamily', 'WindowCount', 'SelectedCoverage', ...
    'AdequateCoverage', 'MedianSelectedN', 'DominantFamily', ...
    'DominantFamilyFraction', 'MedianCondition', 'MedianSigmaMin', ...
    'MedianGeometricVolume', 'MedianEquationResidual', ...
    'MedianPowerResidual', 'MedianEnergyResidual', ...
    'MaxBranchRelativeError', 'Pass'});
selection = vertcat(selectionTables{:});
candidates = vertcat(candidateTables{:});
parameters = vertcat(parameterTables{:});

prefix = 'three_wire_model_order';
if strcmpi(mode, 'quick')
    prefix = 'three_wire_model_order_quick';
end
writetable(summary, fullfile(outDir, [prefix '_summary.csv']));
writetable(selection, fullfile(outDir, [prefix '_selection.csv']));
writetable(candidates, fullfile(outDir, [prefix '_candidates.csv']));
writetable(parameters, fullfile(outDir, [prefix '_parameters.csv']));
save(fullfile(outDir, [prefix '.mat']), 'summary', 'selection', ...
    'candidates', 'parameters', 'config');
plot_model_order_summary(summary, outDir, prefix);

campaign = struct();
campaign.summary = summary;
campaign.selection = selection;
campaign.candidates = candidates;
campaign.parameters = parameters;
campaign.config = config;
campaign.outputPrefix = prefix;

disp(summary);
end

function [vabComponents, vbcComponents] = excitation_components(name)
switch name
    case 'fundamental_balanced'
        vabComponents = [100, 1, 0.00];
        vbcComponents = [100, 1, -2*pi/3];
    case 'fundamental_unbalanced'
        vabComponents = [100, 1, 0.22];
        vbcComponents = [83, 1, -1.86];
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

function truth = truth_for_family(name)
g0 = 1 / 20;
gamma0 = 1 / 0.080;
c0 = 80e-6;
switch name
    case 'balanced'
        g = g0 * [1, 1, 1];
        gamma = gamma0 * [1, 1, 1];
        c = c0 * [1, 1, 1];
    case 'pair_ab_bc'
        g = g0 * [1, 1, 1.45];
        gamma = gamma0 * [1, 1, 0.72];
        c = c0 * [1, 1, 1.35];
    case 'full_unbalanced'
        g = g0 * [1, 1.55, 0.62];
        gamma = gamma0 * [1, 0.64, 1.75];
        c = c0 * [1, 1.42, 0.58];
    otherwise
        error('Unknown truth family "%s".', name);
end
truth = struct();
truth.Gab = g(1);
truth.Gbc = g(2);
truth.Gca = g(3);
truth.Gammaab = gamma(1);
truth.Gammabc = gamma(2);
truth.Gammaca = gamma(3);
truth.Cab = c(1);
truth.Cbc = c(2);
truth.Cca = c(3);
truth.Rab = 1 / truth.Gab;
truth.Rbc = 1 / truth.Gbc;
truth.Rca = 1 / truth.Gca;
truth.Lab = 1 / truth.Gammaab;
truth.Lbc = 1 / truth.Gammabc;
truth.Lca = 1 / truth.Gammaca;
end

function data = make_case(config, truth, vabComponents, vbcComponents, degradationName, caseIndex)
opts = struct();
opts.model = 'delta_parallel_glc';
opts.truth = truth;
opts.vabComponents = vabComponents;
opts.vbcComponents = vbcComponents;
opts.randomSeed = 800 + caseIndex;
switch degradationName
    case 'clean'
        opts.noiseSnrDb = Inf;
    case 'noise_60db'
        opts.noiseSnrDb = 60;
    case 'noise_40db'
        opts.noiseSnrDb = 40;
    otherwise
        error('Unknown degradation "%s".', degradationName);
end
data = nv3_synthetic_delta_parallel_case(config, opts);
data.degradationName = degradationName;
end

function family = expected_family(excitationName, truthName)
if startsWith(excitationName, 'fundamental')
    family = 'reject';
    return;
end
switch truthName
    case 'balanced'
        family = 'balanced';
    case 'pair_ab_bc'
        family = 'pair_equal';
    case 'full_unbalanced'
        family = 'full_unbalanced';
    otherwise
        error('Unknown truth "%s".', truthName);
end
end

function yes = case_passes(result, expectedFamily, dominantFamily, dominantFraction)
if strcmp(expectedFamily, 'reject')
    yes = result.summary.SelectedCoverage <= 0.20;
else
    yes = result.summary.SelectedCoverage >= 0.80 && ...
        result.summary.AdequateCoverage >= 0.80 && ...
        strcmp(dominantFamily, expectedFamily) && ...
        dominantFraction >= 0.80 && ...
        result.summary.MaxBranchRelativeError <= 0.05;
end
end

function family = dominant_selected_family(windows)
selected = windows.Selected > 0;
if ~any(selected)
    family = 'reject';
    return;
end
families = {'balanced', 'pair_equal', 'full_unbalanced'};
counts = zeros(size(families));
for idx = 1:numel(families)
    counts(idx) = sum(strcmp(windows.SelectedFamily(selected), families{idx}));
end
[~, maxIdx] = max(counts);
family = families{maxIdx};
end

function fraction = dominant_selected_fraction(windows, family)
if strcmp(family, 'reject')
    fraction = mean(windows.Selected == 0);
    return;
end
selected = windows.Selected > 0;
if ~any(selected)
    fraction = 0;
else
    fraction = mean(strcmp(windows.SelectedFamily(selected), family));
end
end

function plot_model_order_summary(summary, outDir, prefix)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 900]);
tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
labels = compose('%d %s %s %s', summary.CaseIndex, string(summary.Excitation), ...
    string(summary.TruthFamily), string(summary.Degradation));
x = 1:height(summary);

nexttile;
bar(x, summary.SelectedCoverage);
grid on;
ylim([0 1.05]);
ylabel('selected coverage');
title('Accepted geometric windows');

nexttile;
bar(x, summary.MedianSelectedN);
grid on;
ylim([0 10]);
ylabel('median n_\theta');
title('Selected model order');

nexttile;
semilogy(x, summary.MaxBranchRelativeError, 'o-', 'LineWidth', 1.1);
grid on;
ylabel('max relative error');
title('Branch-parameter error');

nexttile;
passValue = double(summary.Pass);
bar(x, passValue);
grid on;
ylim([0 1.05]);
ylabel('pass');
title('Expected family check');

for ax = findall(fig, 'Type', 'axes').'
    ax.XTick = x;
    ax.XTickLabel = labels;
    ax.XTickLabelRotation = 70;
    ax.FontSize = 8;
end
exportgraphics(fig, fullfile(outDir, [prefix '_summary.png']), 'Resolution', 160);
close(fig);
end
