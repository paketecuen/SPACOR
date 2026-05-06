function campaign = run_three_wire_stress_campaign(mode)
%RUN_THREE_WIRE_STRESS_CAMPAIGN Stress sweep for 3-wire windowed method.

root = fileparts(mfilename('fullpath'));
addpath(root);
config = nv3_default_config();
if ~exist(config.resultsDir, 'dir')
    mkdir(config.resultsDir);
end
if nargin < 1 || isempty(mode)
    mode = 'full';
end

models = {'delta_parallel_gl', 'wye_series_rl'};
[scenarios, noiseLevels, windowCycles, degradations, outputPrefix, makePlots] = ...
    campaign_grid(mode);

rows = {};
allParameters = {};
caseCounter = 0;

for modelIdx = 1:numel(models)
    model = models{modelIdx};
    for scenarioIdx = 1:numel(scenarios)
        scenario = scenarios{scenarioIdx};
        for noiseIdx = 1:numel(noiseLevels)
            noiseSnrDb = noiseLevels(noiseIdx);
            for degradationIdx = 1:numel(degradations)
                degradationName = degradations{degradationIdx};
                for windowIdx = 1:numel(windowCycles)
                    wc = windowCycles(windowIdx);
                    caseCounter = caseCounter + 1;
                    caseName = sprintf('%s_%s_snr%s_%s_w%g', model, scenario, ...
                        snr_label(noiseSnrDb), degradationName, wc);

                    data = make_case(model, scenario, noiseSnrDb, degradationName, config);
                    result = nv3_identify_three_wire_windowed(data, config, struct( ...
                        'harmonics', 'adaptive', ...
                        'windowCycles', wc, ...
                        'hopCycles', min(0.1, wc / 2)));

                    [summaryRow, parameterTable] = summarize_case(caseName, model, ...
                        scenario, noiseSnrDb, degradationName, wc, result);
                    rows(end + 1, :) = summaryRow; %#ok<AGROW>
                    allParameters{end + 1, 1} = parameterTable; %#ok<AGROW>

                    if makePlots && should_plot_case(scenario, noiseSnrDb, degradationName, wc)
                        nv3_plot_parameter_windows(result, data, ...
                            fullfile(config.resultsDir, [caseName '_parameters.png']), ...
                            struct('truthMargin', config.plotTruthMargin));
                    end

                    if mod(caseCounter, 12) == 0
                        fprintf('Completed %d stress cases\n', caseCounter);
                    end
                end
            end
        end
    end
end

summary = cell2table(rows, 'VariableNames', {'Case', 'Model', 'Scenario', ...
    'NoiseSnrDb', 'Degradation', 'WindowCycles', 'ValidCoverage', ...
    'MedianCondition', 'MedianEquationResidual', 'MedianPowerResidual', ...
    'MaxRelativeError', 'PassGate'});
parameters = vertcat(allParameters{:});

writetable(summary, fullfile(config.resultsDir, [outputPrefix '_summary.csv']));
writetable(parameters, fullfile(config.resultsDir, [outputPrefix '_parameters.csv']));
if makePlots
    plot_heatmaps(summary, config, outputPrefix);
end
save(fullfile(config.resultsDir, [outputPrefix '_campaign.mat']), ...
    'summary', 'parameters', 'config');

campaign = struct();
campaign.summary = summary;
campaign.parameters = parameters;
campaign.config = config;

disp(summary);

end

function [scenarios, noiseLevels, windowCycles, degradations, outputPrefix, makePlots] = campaign_grid(mode)
switch lower(char(mode))
    case 'full'
        scenarios = {'fundamental', 'weak_harmonics', 'rich_harmonics'};
        noiseLevels = [Inf, 60, 40];
        windowCycles = [0.25, 0.5, 1.0, 2.0];
        degradations = {'none', 'current_delay_50us'};
        outputPrefix = 'three_wire_stress';
        makePlots = true;
    case 'quick'
        scenarios = {'fundamental', 'rich_harmonics'};
        noiseLevels = [Inf, 60];
        windowCycles = [0.5, 1.0];
        degradations = {'none', 'current_delay_50us'};
        outputPrefix = 'three_wire_stress_quick';
        makePlots = false;
    otherwise
        error('Unknown stress campaign mode "%s". Use "full" or "quick".', char(mode));
end
end

function data = make_case(model, scenario, noiseSnrDb, degradationName, config)
opts = scenario_options(model, scenario);
opts.noiseSnrDb = noiseSnrDb;
switch model
    case 'delta_parallel_gl'
        data = nv3_synthetic_delta_parallel_case(config, opts);
    case 'wye_series_rl'
        data = nv3_synthetic_wye_series_case(config, opts);
    otherwise
        error('Unknown model "%s".', model);
end
data = nv3_apply_measurement_degradation(data, degradation_spec(degradationName));
data.scenario = scenario;
end

function opts = scenario_options(model, scenario)
opts = struct();
switch model
    case 'delta_parallel_gl'
        switch scenario
            case 'fundamental'
                opts.vabComponents = [100, 1, 0.00];
                opts.vbcComponents = [92, 1, -2.18];
            case 'weak_harmonics'
                opts.vabComponents = [100, 1, 0.00; 0.60, 5, 0.40; 0.30, 7, -1.00];
                opts.vbcComponents = [92, 1, -2.18; 0.45, 5, -0.75; 0.25, 7, 1.15];
            case 'rich_harmonics'
                opts.vabComponents = [100, 1, 0.00; 8, 5, 0.40; 4, 7, -1.00];
                opts.vbcComponents = [92, 1, -2.18; 6, 5, -0.75; 3, 7, 1.15];
        end
    case 'wye_series_rl'
        switch scenario
            case 'fundamental'
                opts.iaComponents = [10.0, 1, 0.10];
                opts.ibComponents = [8.0, 1, -2.05];
            case 'weak_harmonics'
                opts.iaComponents = [10.0, 1, 0.10; 0.10, 5, -0.30; 0.05, 7, 1.20];
                opts.ibComponents = [8.0, 1, -2.05; 0.08, 5, 0.85; 0.04, 7, -1.40];
            case 'rich_harmonics'
                opts.iaComponents = [10.0, 1, 0.10; 0.90, 5, -0.30; 0.45, 7, 1.20];
                opts.ibComponents = [8.0, 1, -2.05; 0.70, 5, 0.85; 0.25, 7, -1.40];
        end
end
end

function degradation = degradation_spec(name)
switch char(name)
    case 'none'
        degradation = struct();
    case 'current_delay_50us'
        degradation = struct('currentDelaySec', 50e-6);
    otherwise
        error('Unknown degradation "%s".', char(name));
end
end

function label = snr_label(snrDb)
if isfinite(snrDb)
    label = sprintf('%g', snrDb);
else
    label = 'inf';
end
end

function [row, parameterTable] = summarize_case(caseName, model, scenario, noiseSnrDb, ...
    degradationName, windowCycles, result)
userData = result.summary.Properties.UserData;
finiteError = result.summary.RelativeError(isfinite(result.summary.RelativeError));
if isempty(finiteError)
    maxRelError = NaN;
else
    maxRelError = max(finiteError);
end
passGate = userData.ValidCoverage >= 0.80 && maxRelError <= 0.05 && ...
    userData.MedianPowerResidual <= 0.05 && userData.MedianCondition <= 1e6;
row = {caseName, model, scenario, noiseSnrDb, degradationName, windowCycles, ...
    userData.ValidCoverage, userData.MedianCondition, ...
    userData.MedianEquationResidual, userData.MedianPowerResidual, ...
    maxRelError, passGate};

parameterTable = result.summary;
parameterTable.Case = repmat({caseName}, height(parameterTable), 1);
parameterTable.Model = repmat({model}, height(parameterTable), 1);
parameterTable.Scenario = repmat({scenario}, height(parameterTable), 1);
parameterTable.NoiseSnrDb = repmat(noiseSnrDb, height(parameterTable), 1);
parameterTable.Degradation = repmat({degradationName}, height(parameterTable), 1);
parameterTable.WindowCycles = repmat(windowCycles, height(parameterTable), 1);
parameterTable = movevars(parameterTable, ...
    {'Case', 'Model', 'Scenario', 'NoiseSnrDb', 'Degradation', 'WindowCycles'}, ...
    'Before', 1);
end

function tf = should_plot_case(scenario, noiseSnrDb, degradationName, windowCycles)
tf = strcmp(scenario, 'rich_harmonics') && noiseSnrDb == 60 && ...
    strcmp(degradationName, 'none') && any(abs(windowCycles - [0.25, 0.5, 1.0]) < eps);
end

function plot_heatmaps(summary, config, outputPrefix)
models = unique(summary.Model, 'stable');
scenarios = unique(summary.Scenario, 'stable');
degradations = unique(summary.Degradation, 'stable');
for modelIdx = 1:numel(models)
    model = models{modelIdx};
    for scenarioIdx = 1:numel(scenarios)
        scenario = scenarios{scenarioIdx};
        for degradationIdx = 1:numel(degradations)
            degradation = degradations{degradationIdx};
            mask = strcmp(summary.Model, model) & strcmp(summary.Scenario, scenario) & ...
                strcmp(summary.Degradation, degradation);
            if ~any(mask)
                continue;
            end
            fileBase = sprintf('%s_%s_%s_%s', outputPrefix, model, scenario, degradation);
            plot_metric_heatmap(summary(mask, :), 'MaxRelativeError', ...
                [fileBase '_max_rel_error.png'], config);
            plot_metric_heatmap(summary(mask, :), 'MedianPowerResidual', ...
                [fileBase '_power_residual.png'], config);
            plot_metric_heatmap(summary(mask, :), 'MedianCondition', ...
                [fileBase '_condition.png'], config);
        end
    end
end
end

function plot_metric_heatmap(rows, metricName, fileName, config)
noiseLevels = unique(rows.NoiseSnrDb);
noiseLevels = sort_noise(noiseLevels);
windows = unique(rows.WindowCycles);
windows = sort(windows(:).');
Z = NaN(numel(noiseLevels), numel(windows));
for i = 1:numel(noiseLevels)
    for j = 1:numel(windows)
        mask = rows.NoiseSnrDb == noiseLevels(i) & rows.WindowCycles == windows(j);
        if any(mask)
            Z(i, j) = rows.(metricName)(find(mask, 1));
        end
    end
end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 900, 520]);
imagesc(windows, 1:numel(noiseLevels), log10(max(Z, realmin)));
set(gca, 'YTick', 1:numel(noiseLevels), 'YTickLabel', noise_labels(noiseLevels));
xlabel('Window (cycles)');
ylabel('SNR (dB)');
title(strrep(metricName, '_', ' '));
cb = colorbar;
cb.Label.String = ['log10(' metricName ')'];
grid on;
exportgraphics(fig, fullfile(config.resultsDir, fileName), 'Resolution', 160);
close(fig);
end

function values = sort_noise(values)
finiteValues = sort(values(isfinite(values)), 'descend');
if any(~isfinite(values))
    values = [Inf; finiteValues(:)];
else
    values = finiteValues(:);
end
end

function labels = noise_labels(values)
labels = cell(numel(values), 1);
for idx = 1:numel(values)
    if isfinite(values(idx))
        labels{idx} = sprintf('%g', values(idx));
    else
        labels{idx} = 'Inf';
    end
end
end
