function campaign = run_singlephase_quick(config, opts)
%RUN_SINGLEPHASE_QUICK Canonical quick single-phase SPACOR campaign.
%
% The campaign is intentionally small and reproducible.  It exercises the
% canonical single-phase API on physically transparent synthetic cases and
% writes the standard result bundle under matlab/results/singlephase.

spacor.core.add_legacy_paths();
if nargin < 1 || isempty(config)
    config = spacor.singlephase.default_config();
end
if nargin < 2
    opts = struct();
end

config.fs = spacor.core.get_option(config, 'fs', 10000);
config.f0 = spacor.core.get_option(config, 'f0', 50);
config.duration = spacor.core.get_option(config, 'duration', 0.25);
config.windowCycles = spacor.core.get_option(config, 'windowCycles', 0.5);
config.hopCycles = spacor.core.get_option(config, 'hopCycles', 0.10);
config.harmonics = spacor.core.get_option(config, 'harmonics', [1 3 5 7]);
config.resultsDir = spacor.core.get_option(config, 'resultsDir', ...
    fullfile(spacor.core.lab_root(), 'results', 'singlephase'));

outDir = spacor.core.get_option(opts, 'outputDir', ...
    fullfile(config.resultsDir, 'quick'));
outDir = spacor.io.ensure_dir(outDir);
figureDir = spacor.io.ensure_dir(fullfile(outDir, 'figures'));
candidateMode = spacor.core.get_option(opts, 'candidateMode', 'lti');

caseDefs = singlephase_cases();
summaryRows = cell(numel(caseDefs), 11);
parameterTables = cell(numel(caseDefs), 1);
stabilityRows = cell(numel(caseDefs), 8);
informationRows = cell(numel(caseDefs), 8);
rankingTables = cell(numel(caseDefs), 1);
reports = cell(numel(caseDefs), 1);
dataCases = cell(numel(caseDefs), 1);

for idx = 1:numel(caseDefs)
    spec = caseDefs(idx);
    data = make_singlephase_case(spec, config);
    report = spacor.singlephase.characterize(data, config, ...
        struct('candidateMode', candidateMode));
    reports{idx} = report;
    dataCases{idx} = data;
    spacor.plot.parameter_evolution(report, data, ...
        fullfile(figureDir, [safe_name(spec.name) '_parameters.png']), ...
        struct('truthMargin', config.plotTruthMargin));

    selectedIdx = selected_candidate_index(report);
    [residual, coverage, confidence] = selected_metrics(report, selectedIdx);
    pass = is_pass(report, coverage);

    summaryRows(idx, :) = {idx, spec.name, spec.truthModel, ...
        report.canonical.status, report.canonical.selectedModel, ...
        report.canonical.family, pass, residual, coverage, confidence, ...
        report.recommendation.reason};

    parameterTables{idx} = parameter_table(idx, spec, report, selectedIdx);
    stabilityRows(idx, :) = {idx, spec.name, pass, residual, coverage, ...
        confidence, report.input.vRms, report.input.iRms};
    informationRows(idx, :) = {idx, spec.name, data.fs, data.f0, ...
        numel(data.t), data.t(end) - data.t(1), ...
        rms_value(data.v), rms_value(data.i)};

    ranking = report.candidates;
    ranking.CaseIndex = repmat(idx, height(ranking), 1);
    ranking.CaseName = repmat({spec.name}, height(ranking), 1);
    ranking = movevars(ranking, {'CaseIndex', 'CaseName'}, 'Before', 1);
    rankingTables{idx} = ranking;
end

summary = cell2table(summaryRows, 'VariableNames', { ...
    'CaseIndex', 'CaseName', 'TruthModel', 'Status', 'SelectedModel', ...
    'SelectedFamily', 'Pass', 'ResidualOrScore', 'Coverage', ...
    'Confidence', 'Reason'});
parameters = vertcat(parameterTables{:});
stability = cell2table(stabilityRows, 'VariableNames', { ...
    'CaseIndex', 'CaseName', 'Pass', 'ResidualOrScore', 'Coverage', ...
    'Confidence', 'Vrms', 'Irms'});
information = cell2table(informationRows, 'VariableNames', { ...
    'CaseIndex', 'CaseName', 'Fs', 'F0', 'Samples', 'Duration', ...
    'Vrms', 'Irms'});
method_ranking = vertcat(rankingTables{:});
summaryFigure = fullfile(figureDir, 'campaign_summary.png');
spacor.plot.campaign_summary(struct('summary', summary, ...
    'stability', stability, 'outputDir', outDir), summaryFigure);

writetable(summary, fullfile(outDir, 'summary.csv'));
writetable(parameters, fullfile(outDir, 'parameters.csv'));
writetable(stability, fullfile(outDir, 'stability.csv'));
writetable(information, fullfile(outDir, 'information.csv'));
writetable(method_ranking, fullfile(outDir, 'method_ranking.csv'));
save(fullfile(outDir, 'report.mat'), 'summary', 'parameters', ...
    'stability', 'information', 'method_ranking', 'reports', ...
    'dataCases', 'config', 'opts', 'summaryFigure');

campaign = struct('summary', summary, 'parameters', parameters, ...
    'stability', stability, 'information', information, ...
    'method_ranking', method_ranking, 'reports', {reports}, ...
    'dataCases', {dataCases}, 'config', config, 'outputDir', outDir, ...
    'figureDir', figureDir, 'summaryFigure', summaryFigure, ...
    'api', struct('name', 'spacor.campaigns.run_singlephase_quick'));
end

function cases = singlephase_cases()
cases = struct('name', {}, 'truthModel', {}, 'kind', {});
cases(end + 1) = struct('name', 'series_rl_multisine', ...
    'truthModel', 'lti_series_rl', 'kind', 'series_rl');
cases(end + 1) = struct('name', 'parallel_rc_multisine', ...
    'truthModel', 'lti_parallel_rc', 'kind', 'parallel_rc');
cases(end + 1) = struct('name', 'parallel_series_rl_multisine', ...
    'truthModel', 'lti_parallel_series_rl', 'kind', 'parallel_series_rl');
end

function data = make_singlephase_case(spec, config)
t = (0:round(config.duration * config.fs) - 1).' / config.fs;
w0 = 2 * pi * config.f0;
v = multisine(t, w0, [120 9 5 3], [1 3 5 7], [0 0.4 -0.7 1.1]);

switch spec.kind
    case 'series_rl'
        truth = struct('R', 10, 'L', 50e-3);
        i = current_series_rl(t, w0, truth, [120 9 5 3], ...
            [1 3 5 7], [0 0.4 -0.7 1.1]);
    case 'parallel_rc'
        truth = struct('G', 1 / 25, 'R', 25, 'C', 120e-6);
        dv = multisine_derivative(t, w0, [120 9 5 3], ...
            [1 3 5 7], [0 0.4 -0.7 1.1], 1);
        i = truth.G * v + truth.C * dv;
    case 'parallel_series_rl'
        truth = struct('Gp', 1 / 8, 'Rp', 8, 'Rser', 2.5, ...
            'Lser', 70e-3);
        iSeries = current_series_rl(t, w0, ...
            struct('R', truth.Rser, 'L', truth.Lser), ...
            [120 9 5 3], [1 3 5 7], [0 0.4 -0.7 1.1]);
        i = truth.Gp * v + iSeries;
    otherwise
        error('Unknown single-phase quick case "%s".', spec.kind);
end

data = struct('t', t, 'v', v, 'i', i, 'fs', config.fs, ...
    'f0', config.f0, 'label', spec.name, 'truth', truth, ...
    'truthModel', spec.truthModel);
end

function y = multisine(t, w0, amps, harmonics, phases)
y = zeros(size(t));
for idx = 1:numel(harmonics)
    omega = harmonics(idx) * w0;
    y = y + amps(idx) * sin(omega * t + phases(idx));
end
end

function y = multisine_derivative(t, w0, amps, harmonics, phases, order)
y = zeros(size(t));
for idx = 1:numel(harmonics)
    omega = harmonics(idx) * w0;
    y = y + amps(idx) * omega^order * ...
        sin(omega * t + phases(idx) + order * pi / 2);
end
end

function i = current_series_rl(t, w0, truth, amps, harmonics, phases)
i = zeros(size(t));
for idx = 1:numel(harmonics)
    omega = harmonics(idx) * w0;
    z = truth.R + 1i * omega * truth.L;
    i = i + amps(idx) / abs(z) * ...
        sin(omega * t + phases(idx) - angle(z));
end
end

function idx = selected_candidate_index(report)
idx = find(strcmp(report.candidates.Candidate, report.canonical.selectedModel), 1);
if isempty(idx) && strcmp(report.canonical.selectedModel, 'lti_equivalent_family')
    accepted = report.candidates.Accepted > 0 & ...
        strcmp(report.candidates.Family, 'lti_windowed');
    if any(accepted)
        [~, localIdx] = min(report.candidates.ResidualOrScore(accepted));
        indices = find(accepted);
        idx = indices(localIdx);
    end
end
end

function [residual, coverage, confidence] = selected_metrics(report, idx)
if isempty(idx)
    residual = NaN;
    coverage = NaN;
    confidence = NaN;
else
    residual = report.candidates.ResidualOrScore(idx);
    coverage = report.candidates.Coverage(idx);
    confidence = report.candidates.Confidence(idx);
end
end

function pass = is_pass(report, coverage)
pass = startsWith(string(report.canonical.status), "accepted") && ...
    isfinite(coverage) && coverage >= 0.5;
end

function tableOut = parameter_table(caseIndex, spec, report, selectedIdx)
if isempty(selectedIdx)
    params = "";
else
    params = string(report.candidates.Parameters{selectedIdx});
end
tableOut = table(caseIndex, string(spec.name), string(spec.truthModel), ...
    string(report.canonical.selectedModel), params, ...
    'VariableNames', {'CaseIndex', 'CaseName', 'TruthModel', ...
    'SelectedModel', 'ParameterSummary'});
end

function value = rms_value(x)
value = sqrt(mean(x(:) .^ 2));
end

function textOut = safe_name(textIn)
textOut = regexprep(char(textIn), '[^A-Za-z0-9_]+', '_');
textOut = regexprep(textOut, '_+', '_');
textOut = regexprep(textOut, '^_|_$', '');
end
