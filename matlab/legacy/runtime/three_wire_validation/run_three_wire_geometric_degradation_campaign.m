function campaign = run_three_wire_geometric_degradation_campaign(mode)
%RUN_THREE_WIRE_GEOMETRIC_DEGRADATION_CAMPAIGN Degrade full delta RLC tests.
%
% This campaign uses the practical windowed estimator with fitted
% time-domain primitives and derivatives. It tests the fully unbalanced
% delta parallel RLC model:
%
%   i_xy = G_xy v_xy + Gamma_xy integral(v_xy) + C_xy dv_xy/dt.

if nargin < 1 || isempty(mode)
    mode = 'full';
end

config = nv3_default_config();
config.duration = 0.8;
config.fs = 20000;
config.f0 = 50;
config.adaptiveHarmonics = true;
config.harmonicCandidates = [1 3 5 7 9 11 13];
config.activeThreshold = 1e-3;
config.conditionLimit = 1e8;

if strcmpi(mode, 'quick')
    windowCyclesList = [0.5, 1.0];
    degradationCases = base_degradation_cases();
else
    windowCyclesList = [0.25, 0.5, 1.0, 2.0];
    degradationCases = [base_degradation_cases(), extra_degradation_cases()];
end

outDir = config.resultsDir;
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

baseOpts = struct();
baseOpts.model = 'delta_parallel_glc';
baseOpts.noiseSnrDb = Inf;
baseOpts.randomSeed = 11;
baseOpts.vabComponents = [100, 1, 0.00; 8, 5, 0.40; 4, 7, -1.00];
baseOpts.vbcComponents = [92, 1, -2.18; 6, 5, -0.75; 3, 7, 1.15];

summaryRows = {};
parameterTables = {};
caseIndex = 0;
representativeDone = false;

for d = 1:numel(degradationCases)
    degCase = degradationCases(d);
    for wc = 1:numel(windowCyclesList)
        caseIndex = caseIndex + 1;
        opts = baseOpts;
        opts.noiseSnrDb = degCase.noiseSnrDb;
        opts.randomSeed = 101 + d;
        data = nv3_synthetic_delta_parallel_case(config, opts);
        if isfield(degCase, 'degradation')
            data = nv3_apply_measurement_degradation(data, degCase.degradation);
        end
        if isfield(degCase, 'quantizationBits') && isfinite(degCase.quantizationBits)
            data = quantize_measurements(data, degCase.quantizationBits);
        end
        idOpts = struct('model', 'delta_parallel_glc', ...
            'windowCycles', windowCyclesList(wc), ...
            'hopCycles', 0.1, ...
            'adaptiveHarmonics', true);
        result = nv3_identify_three_wire_windowed(data, config, idOpts);
        [maxErr, p90Err] = error_stats(result.summary);
        validCoverage = result.summary.Properties.UserData.ValidCoverage;
        medianCondition = result.summary.Properties.UserData.MedianCondition;
        medianPowerResidual = result.summary.Properties.UserData.MedianPowerResidual;
        medianEquationResidual = result.summary.Properties.UserData.MedianEquationResidual;
        pass = validCoverage >= 0.80 && maxErr <= 0.05 && ...
            p90Err <= 0.15 && medianCondition <= config.conditionLimit && ...
            medianPowerResidual <= 0.05;
        summaryRows(end + 1, :) = {caseIndex, degCase.name, ...
            windowCyclesList(wc), degCase.noiseSnrDb, validCoverage, ...
            medianCondition, medianEquationResidual, medianPowerResidual, ...
            maxErr, p90Err, pass}; %#ok<AGROW>
        params = result.summary;
        params.CaseIndex = repmat(caseIndex, height(params), 1);
        params.Degradation = repmat(string(degCase.name), height(params), 1);
        params.WindowCycles = repmat(windowCyclesList(wc), height(params), 1);
        parameterTables{end + 1} = params; %#ok<AGROW>

        if should_plot_representative(degCase.name, windowCyclesList(wc), representativeDone)
            plotName = sprintf('three_wire_delta_rlc_%s_w%.2g_waveforms_parameters.png', ...
                degCase.name, windowCyclesList(wc));
            nv3_plot_waveforms_and_parameters(result, data, fullfile(outDir, plotName), ...
                struct('truthMargin', 0.08));
            representativeDone = representativeDone || strcmp(degCase.name, 'clean');
        end
    end
end

summary = cell2table(summaryRows, 'VariableNames', { ...
    'CaseIndex', 'Degradation', 'WindowCycles', 'NoiseSnrDb', ...
    'Coverage', 'MedianCondition', 'MedianEquationResidual', ...
    'MedianPowerResidual', 'MaxRelativeError', 'P90RelativeError', 'Pass'});
parameters = vertcat(parameterTables{:});

prefix = 'three_wire_geometric_degradation';
if strcmpi(mode, 'quick')
    prefix = 'three_wire_geometric_degradation_quick';
end
writetable(summary, fullfile(outDir, [prefix '_summary.csv']));
writetable(parameters, fullfile(outDir, [prefix '_parameters.csv']));
save(fullfile(outDir, [prefix '.mat']), 'summary', 'parameters');
plot_degradation_summary(summary, outDir, prefix);

campaign = struct('summary', summary, 'parameters', parameters, ...
    'config', config, 'outputPrefix', prefix);
disp(summary);
end

function cases = base_degradation_cases()
cases = struct('name', {}, 'noiseSnrDb', {}, 'degradation', {}, 'quantizationBits', {});
cases(end + 1) = struct('name', 'clean', 'noiseSnrDb', Inf, ...
    'degradation', struct(), 'quantizationBits', Inf);
cases(end + 1) = struct('name', 'noise_60db', 'noiseSnrDb', 60, ...
    'degradation', struct(), 'quantizationBits', Inf);
cases(end + 1) = struct('name', 'noise_40db', 'noiseSnrDb', 40, ...
    'degradation', struct(), 'quantizationBits', Inf);
cases(end + 1) = struct('name', 'current_delay_50us', 'noiseSnrDb', Inf, ...
    'degradation', struct('currentDelaySec', 50e-6), 'quantizationBits', Inf);
end

function cases = extra_degradation_cases()
cases = struct('name', {}, 'noiseSnrDb', {}, 'degradation', {}, 'quantizationBits', {});
cases(end + 1) = struct('name', 'gain_mismatch', 'noiseSnrDb', Inf, ...
    'degradation', struct('voltageGain', [1.002 0.998 1.000], ...
    'currentGain', [0.997 1.003 1.000]), 'quantizationBits', Inf);
cases(end + 1) = struct('name', 'delay_plus_noise_60db', 'noiseSnrDb', 60, ...
    'degradation', struct('currentDelaySec', 50e-6), 'quantizationBits', Inf);
cases(end + 1) = struct('name', 'quantized_12bit', 'noiseSnrDb', Inf, ...
    'degradation', struct(), 'quantizationBits', 12);
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

function [maxErr, p90Err] = error_stats(summary)
mask = isfinite(summary.RelativeError) & summary.Truth > 0;
err = summary.RelativeError(mask);
if isempty(err)
    maxErr = NaN;
    p90Err = NaN;
else
    maxErr = max(err);
    p90Err = percentile(err, 90);
end
end

function yes = should_plot_representative(name, windowCycles, representativeDone)
yes = false;
if abs(windowCycles - 0.5) > 1e-12
    return;
end
if strcmp(name, 'clean') && ~representativeDone
    yes = true;
elseif any(strcmp(name, {'noise_60db', 'current_delay_50us', 'delay_plus_noise_60db'}))
    yes = true;
end
end

function plot_degradation_summary(summary, outDir, prefix)
degradations = unique(summary.Degradation, 'stable');
windows = unique(summary.WindowCycles, 'stable');
err = nan(numel(degradations), numel(windows));
coverage = nan(numel(degradations), numel(windows));
for d = 1:numel(degradations)
    for w = 1:numel(windows)
        mask = strcmp(summary.Degradation, degradations{d}) & summary.WindowCycles == windows(w);
        if any(mask)
            err(d, w) = summary.MaxRelativeError(mask);
            coverage(d, w) = summary.Coverage(mask);
        end
    end
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 480]);
tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
imagesc(coverage, [0 1]);
colorbar;
title('Coverage');
set(gca, 'XTick', 1:numel(windows), 'XTickLabel', compose('%.2g', windows), ...
    'YTick', 1:numel(degradations), 'YTickLabel', degradations, ...
    'TickLabelInterpreter', 'none');
xlabel('Window cycles');
nexttile;
imagesc(log10(err));
colorbar;
title('log10 max relative error');
set(gca, 'XTick', 1:numel(windows), 'XTickLabel', compose('%.2g', windows), ...
    'YTick', 1:numel(degradations), 'YTickLabel', degradations, ...
    'TickLabelInterpreter', 'none');
xlabel('Window cycles');
exportgraphics(fig, fullfile(outDir, [prefix '_heatmap.png']), 'Resolution', 160);
close(fig);
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
