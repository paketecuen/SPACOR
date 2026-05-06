function fig = campaign_summary(campaignOrDir, outFile, opts)
%CAMPAIGN_SUMMARY Plot an aggregate dashboard for a SPACOR campaign.
%
% Input can be either a campaign struct with fields summary/stability or a
% directory containing summary.csv and, optionally, stability.csv.

if nargin < 2
    outFile = '';
end
if nargin < 3
    opts = struct();
end

[summary, stability, label] = load_campaign_tables(campaignOrDir);
maxRows = spacor.core.get_option(opts, 'maxRows', 40);
if height(summary) > maxRows
    summary = summary(1:maxRows, :);
end
if ~isempty(stability) && height(stability) > maxRows
    stability = stability(1:maxRows, :);
end

caseLabels = row_labels(summary);
passValues = column_or_nan(summary, 'Pass');
selectedModels = string_column(summary, 'SelectedModel', "unknown");
residual = first_available_numeric(summary, stability, ...
    {'ResidualOrScore', 'MedianPowerResidual', 'BestUsablePowerResidual'});
coverage = first_available_numeric(summary, stability, ...
    {'Coverage', 'SelectedCoverage', 'AdequateCoverage'});
condition = first_available_numeric(summary, stability, ...
    {'MedianCondition', 'ConditionNumber'});

fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [90, 90, 1550, 1100]);
layout = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
title(layout, sprintf('SPACOR campaign summary | %s', label), ...
    'Interpreter', 'none');

nexttile(layout);
plot_pass_panel(caseLabels, passValues);

nexttile(layout);
plot_model_panel(caseLabels, selectedModels);

nexttile(layout);
plot_residual_condition_panel(caseLabels, residual, condition);

nexttile(layout);
plot_coverage_panel(caseLabels, coverage);

if ~isempty(outFile)
    outDir = fileparts(outFile);
    if ~isempty(outDir)
        spacor.io.ensure_dir(outDir);
    end
    exportgraphics(fig, outFile, 'Resolution', 170);
end
end

function [summary, stability, label] = load_campaign_tables(campaignOrDir)
stability = table();
if isstruct(campaignOrDir)
    summary = campaignOrDir.summary;
    if isfield(campaignOrDir, 'stability')
        stability = campaignOrDir.stability;
    end
    if isfield(campaignOrDir, 'outputDir')
        label = campaignOrDir.outputDir;
    else
        label = 'in-memory campaign';
    end
else
    campaignDir = char(campaignOrDir);
    summary = readtable(fullfile(campaignDir, 'summary.csv'));
    stabilityFile = fullfile(campaignDir, 'stability.csv');
    if isfile(stabilityFile)
        stability = readtable(stabilityFile);
    end
    label = campaignDir;
end
end

function labels = row_labels(summary)
caseName = string_column(summary, 'CaseName', "");
if all(caseName == "")
    caseName = "case_" + string((1:height(summary)).');
end
if ismember('Degradation', summary.Properties.VariableNames)
    deg = string_column(summary, 'Degradation', "");
    labels = caseName + " | " + deg;
else
    labels = caseName;
end
labels = matlab.lang.makeUniqueStrings(cellstr(labels));
labels = string(labels(:));
end

function plot_pass_panel(labels, passValues)
x = 1:numel(labels);
if all(~isfinite(passValues))
    text(0.5, 0.5, 'No Pass column available', ...
        'HorizontalAlignment', 'center');
    axis off;
    return;
end
passValues = double(passValues(:) > 0);
b = bar(x, passValues, 'FaceColor', 'flat');
b.CData = repmat([0.80 0.20 0.20], numel(x), 1);
b.CData(passValues > 0, :) = repmat([0.20 0.60 0.30], ...
    sum(passValues > 0), 1);
ylim([0 1.15]);
grid on;
set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 35, ...
    'TickLabelInterpreter', 'none');
ylabel('pass');
title(sprintf('Pass rate %.1f%%', 100 * mean(passValues > 0, 'omitnan')));
end

function plot_model_panel(labels, models)
x = 1:numel(labels);
[uniqueModels, ~, idx] = unique(models, 'stable');
plot(x, idx, 'o-', 'LineWidth', 1.0, 'MarkerFaceColor', [0.0 0.45 0.74]);
grid on;
ylim([0.5, max(1.5, numel(uniqueModels) + 0.5)]);
set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 35, ...
    'YTick', 1:numel(uniqueModels), 'YTickLabel', uniqueModels, ...
    'TickLabelInterpreter', 'none');
ylabel('selected');
title('Selected physical equivalent');
end

function plot_residual_condition_panel(labels, residual, condition)
x = 1:numel(labels);
hasResidual = any(isfinite(residual));
hasCondition = any(isfinite(condition));
if ~hasResidual && ~hasCondition
    text(0.5, 0.5, 'No residual or condition data available', ...
        'HorizontalAlignment', 'center');
    axis off;
    return;
end
if hasResidual
    yyaxis left;
    semilogy(x, positive_or_nan(residual), 'o-', 'LineWidth', 1.0, ...
        'MarkerFaceColor', [0.0 0.45 0.74]);
    ylabel('residual');
end
if hasCondition
    yyaxis right;
    semilogy(x, positive_or_nan(condition), 's-', 'LineWidth', 1.0, ...
        'MarkerFaceColor', [0.85 0.33 0.10]);
    ylabel('condition');
end
grid on;
set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 35, ...
    'TickLabelInterpreter', 'none');
title('Residual and conditioning');
end

function plot_coverage_panel(labels, coverage)
x = 1:numel(labels);
if all(~isfinite(coverage))
    text(0.5, 0.5, 'No coverage data available', ...
        'HorizontalAlignment', 'center');
    axis off;
    return;
end
bar(x, coverage, 'FaceColor', [0.25 0.50 0.75]);
hold on;
yline(0.80, '--', '80% gate', 'LineWidth', 0.9);
hold off;
ylim([0, 1.05 * max([1; coverage(isfinite(coverage))])]);
grid on;
set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 35, ...
    'TickLabelInterpreter', 'none');
ylabel('coverage');
title('Valid or selected window coverage');
end

function values = first_available_numeric(summary, stability, names)
values = NaN(height(summary), 1);
for idx = 1:numel(names)
    if ismember(names{idx}, summary.Properties.VariableNames)
        values = numeric_column(summary, names{idx});
        return;
    end
end
if ~isempty(stability)
    for idx = 1:numel(names)
        if ismember(names{idx}, stability.Properties.VariableNames)
            values = numeric_column(stability, names{idx});
            return;
        end
    end
end
end

function values = column_or_nan(tableIn, name)
if ismember(name, tableIn.Properties.VariableNames)
    values = numeric_column(tableIn, name);
else
    values = NaN(height(tableIn), 1);
end
end

function values = numeric_column(tableIn, name)
raw = tableIn.(name);
if isnumeric(raw) || islogical(raw)
    values = double(raw(:));
else
    values = str2double(string(raw(:)));
end
end

function values = string_column(tableIn, name, defaultValue)
if ismember(name, tableIn.Properties.VariableNames)
    values = string(tableIn.(name));
else
    values = repmat(defaultValue, height(tableIn), 1);
end
values = values(:);
end

function y = positive_or_nan(x)
y = x(:);
y(~isfinite(y) | y <= 0) = NaN;
end

