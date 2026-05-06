function campaign = run_singlephase_degradation(config, opts)
%RUN_SINGLEPHASE_DEGRADATION Canonical single-phase degradation campaign.

spacor.core.add_legacy_paths();
if nargin < 1 || isempty(config)
    config = spacor.singlephase.default_config();
end
if nargin < 2
    opts = struct();
end

mode = char(spacor.core.get_option(opts, 'mode', 'quick'));
config.fs = spacor.core.get_option(config, 'fs', 10000);
config.f0 = spacor.core.get_option(config, 'f0', 50);
config.duration = spacor.core.get_option(config, 'duration', 0.25);
config.windowCycles = spacor.core.get_option(config, 'windowCycles', 0.5);
config.hopCycles = spacor.core.get_option(config, 'hopCycles', 0.10);
config.harmonics = spacor.core.get_option(config, 'harmonics', [1 3 5 7]);
config.resultsDir = spacor.core.get_option(config, 'resultsDir', ...
    fullfile(spacor.core.lab_root(), 'results', 'singlephase'));

outDir = spacor.core.get_option(opts, 'outputDir', ...
    fullfile(config.resultsDir, ['degradation_' mode]));
outDir = spacor.io.ensure_dir(outDir);
figureDir = spacor.io.ensure_dir(fullfile(outDir, 'figures'));
makeFigures = spacor.core.get_option(opts, 'makeFigures', true);

cases = singlephase_cases(mode);
degradations = degradation_specs(mode);
nRows = numel(cases) * numel(degradations);
summaryRows = cell(nRows, 16);
parameterTables = cell(nRows, 1);
stabilityRows = cell(nRows, 12);
informationRows = cell(nRows, 12);
rankingTables = cell(nRows, 1);
reports = cell(nRows, 1);
dataCases = cell(nRows, 1);

row = 0;
for c = 1:numel(cases)
    baseData = make_singlephase_case(cases(c), config);
    for d = 1:numel(degradations)
        row = row + 1;
        deg = degradations(d);
        rng(8300 + 100 * c + d);
        data = apply_singlephase_degradation(baseData, deg);
        data.label = sprintf('%s | %s', cases(c).name, deg.name);
        report = spacor.singlephase.characterize(data, config, ...
            struct('candidateMode', 'lti'));
        reports{row} = report;
        dataCases{row} = data;

        selectedIdx = selected_candidate_index(report);
        [residual, coverage, confidence] = selected_metrics(report, selectedIdx);
        pass = is_singlephase_pass(report, cases(c).truthModel, coverage);
        summaryRows(row, :) = {row, cases(c).name, deg.name, ...
            cases(c).truthModel, report.canonical.status, ...
            report.canonical.selectedModel, report.canonical.family, pass, ...
            residual, coverage, confidence, deg.noiseSnrDb, ...
            deg.gainError, deg.offsetFraction, deg.quantizationBits, ...
            deg.movingAverageSamples};

        parameterTables{row} = singlephase_parameter_table(row, cases(c), ...
            deg, report, selectedIdx);
        stabilityRows(row, :) = {row, cases(c).name, deg.name, pass, ...
            residual, coverage, confidence, report.input.vRms, ...
            report.input.iRms, deg.noiseSnrDb, deg.offsetFraction, ...
            deg.movingAverageSamples};
        informationRows(row, :) = {row, cases(c).name, deg.name, ...
            data.fs, data.f0, numel(data.t), data.t(end) - data.t(1), ...
            rms_value(data.v), rms_value(data.i), deg.noiseSnrDb, ...
            deg.quantizationBits, deg.gainError};

        ranking = report.candidates;
        ranking.CaseIndex = repmat(row, height(ranking), 1);
        ranking.CaseName = repmat({cases(c).name}, height(ranking), 1);
        ranking.Degradation = repmat({deg.name}, height(ranking), 1);
        ranking = movevars(ranking, {'CaseIndex', 'CaseName', ...
            'Degradation'}, 'Before', 1);
        rankingTables{row} = ranking;

        if makeFigures
            safe = sprintf('%03d_%s_%s_parameters.png', row, ...
                safe_name(cases(c).name), safe_name(deg.name));
            try
                spacor.plot.parameter_evolution(report, data, ...
                    fullfile(figureDir, safe), ...
                    struct('truthMargin', config.plotTruthMargin));
            catch err
                warning('spacor:plotFailed', ...
                    'Single-phase degradation plot failed for %s/%s: %s', ...
                    cases(c).name, deg.name, err.message);
            end
        end
    end
end

summary = cell2table(summaryRows, 'VariableNames', { ...
    'CaseIndex', 'CaseName', 'Degradation', 'TruthModel', 'Status', ...
    'SelectedModel', 'SelectedFamily', 'Pass', 'ResidualOrScore', ...
    'Coverage', 'Confidence', 'NoiseSnrDb', 'GainError', ...
    'OffsetFraction', 'QuantizationBits', 'MovingAverageSamples'});
parameters = vertcat(parameterTables{:});
stability = cell2table(stabilityRows, 'VariableNames', { ...
    'CaseIndex', 'CaseName', 'Degradation', 'Pass', ...
    'ResidualOrScore', 'Coverage', 'Confidence', 'Vrms', 'Irms', ...
    'NoiseSnrDb', 'OffsetFraction', 'MovingAverageSamples'});
information = cell2table(informationRows, 'VariableNames', { ...
    'CaseIndex', 'CaseName', 'Degradation', 'Fs', 'F0', 'Samples', ...
    'Duration', 'Vrms', 'Irms', 'NoiseSnrDb', 'QuantizationBits', ...
    'GainError'});
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
    'dataCases', 'config', 'opts', 'mode', 'summaryFigure');

campaign = struct('summary', summary, 'parameters', parameters, ...
    'stability', stability, 'information', information, ...
    'method_ranking', method_ranking, 'reports', {reports}, ...
    'dataCases', {dataCases}, 'config', config, 'mode', mode, ...
    'outputDir', outDir, 'figureDir', figureDir, ...
    'summaryFigure', summaryFigure, ...
    'api', struct('name', 'spacor.campaigns.run_singlephase_degradation'));
end

function cases = singlephase_cases(mode)
cases = struct('name', {}, 'truthModel', {}, 'kind', {});
cases(end + 1) = struct('name', 'series_rl_multisine', ...
    'truthModel', 'lti_series_rl', 'kind', 'series_rl');
cases(end + 1) = struct('name', 'parallel_series_rl_multisine', ...
    'truthModel', 'lti_parallel_series_rl', 'kind', 'parallel_series_rl');
cases(end + 1) = struct('name', 'parallel_rc_multisine', ...
    'truthModel', 'lti_parallel_rc', 'kind', 'parallel_rc');
if strcmp(mode, 'smoke')
    cases = cases(1);
end
end

function specs = degradation_specs(mode)
specs = struct('name', {}, 'noiseSnrDb', {}, 'gainError', {}, ...
    'offsetFraction', {}, 'quantizationBits', {}, ...
    'movingAverageSamples', {});
specs(end + 1) = make_deg('clean', Inf, 0, 0, 0, 0);
specs(end + 1) = make_deg('noise_60db', 60, 0, 0, 0, 0);
if ~strcmp(mode, 'smoke')
    specs(end + 1) = make_deg('noise_40db', 40, 0, 0, 0, 0);
    specs(end + 1) = make_deg('gain_mismatch_1pct', Inf, 0.01, 0, 0, 0);
    specs(end + 1) = make_deg('offset_1pct', Inf, 0, 0.01, 0, 0);
    specs(end + 1) = make_deg('quantized_12bit', Inf, 0, 0, 12, 0);
    specs(end + 1) = make_deg('moving_average_5', Inf, 0, 0, 0, 5);
end
end

function deg = make_deg(name, noiseSnrDb, gainError, offsetFraction, ...
    quantizationBits, movingAverageSamples)
deg = struct('name', name, 'noiseSnrDb', noiseSnrDb, ...
    'gainError', gainError, 'offsetFraction', offsetFraction, ...
    'quantizationBits', quantizationBits, ...
    'movingAverageSamples', movingAverageSamples);
end

function data = make_singlephase_case(spec, config)
t = (0:round(config.duration * config.fs) - 1).' / config.fs;
w0 = 2 * pi * config.f0;
amps = [120 9 5 3];
harmonics = [1 3 5 7];
phases = [0 0.4 -0.7 1.1];
v = multisine(t, w0, amps, harmonics, phases, 0);

switch spec.kind
    case 'series_rl'
        truth = struct('R', 10, 'L', 50e-3);
        i = current_series_rl(t, w0, truth, amps, harmonics, phases);
    case 'parallel_rc'
        truth = struct('G', 1 / 25, 'R', 25, 'C', 120e-6);
        dv = multisine(t, w0, amps, harmonics, phases, 1);
        i = truth.G * v + truth.C * dv;
    case 'parallel_series_rl'
        truth = struct('Gp', 1 / 8, 'Rp', 8, 'Rser', 2.5, ...
            'Rs', 2.5, 'Lser', 70e-3, 'L', 70e-3);
        iSeries = current_series_rl(t, w0, ...
            struct('R', truth.Rser, 'L', truth.Lser), amps, ...
            harmonics, phases);
        i = truth.Gp * v + iSeries;
    otherwise
        error('Unknown single-phase degradation case "%s".', spec.kind);
end

data = struct('t', t, 'v', v, 'i', i, 'fs', config.fs, ...
    'f0', config.f0, 'label', spec.name, 'truth', truth, ...
    'truthModel', spec.truthModel);
end

function data = apply_singlephase_degradation(data, deg)
v = data.v(:);
i = data.i(:);
if isfinite(deg.noiseSnrDb)
    v = add_noise(v, deg.noiseSnrDb);
    i = add_noise(i, deg.noiseSnrDb);
end
if deg.gainError ~= 0
    v = (1 + deg.gainError) * v;
    i = (1 - deg.gainError) * i;
end
if deg.offsetFraction > 0
    v = v + deg.offsetFraction * rms_value(v);
    i = i + deg.offsetFraction * rms_value(i);
end
if deg.quantizationBits > 0
    v = quantize_signal(v, deg.quantizationBits);
    i = quantize_signal(i, deg.quantizationBits);
end
if deg.movingAverageSamples > 1
    v = movmean(v, deg.movingAverageSamples, 'Endpoints', 'shrink');
    i = movmean(i, deg.movingAverageSamples, 'Endpoints', 'shrink');
end
data.v = v;
data.i = i;
data.degradation = deg;
data.noiseSnrDb = deg.noiseSnrDb;
end

function y = multisine(t, w0, amps, harmonics, phases, order)
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
        indices = find(accepted);
        [~, local] = min(report.candidates.ResidualOrScore(indices));
        idx = indices(local);
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

function pass = is_singlephase_pass(report, truthModel, coverage)
selected = string(report.canonical.selectedModel);
truth = string(truthModel);
accepted = startsWith(string(report.canonical.status), "accepted");
compatible = selected == truth || selected == "lti_equivalent_family";
pass = accepted && compatible && isfinite(coverage) && coverage >= 0.5;
end

function tableOut = singlephase_parameter_table(caseIndex, spec, deg, report, idx)
if isempty(idx)
    params = "";
else
    params = string(report.candidates.Parameters{idx});
end
tableOut = table(caseIndex, string(spec.name), string(deg.name), ...
    string(spec.truthModel), string(report.canonical.selectedModel), params, ...
    'VariableNames', {'CaseIndex', 'CaseName', 'Degradation', ...
    'TruthModel', 'SelectedModel', 'ParameterSummary'});
end

function y = add_noise(x, snrDb)
sigma = rms_value(x) * 10 ^ (-snrDb / 20);
y = x + sigma * randn(size(x));
end

function y = quantize_signal(x, bits)
levels = 2 ^ bits;
span = max(x) - min(x);
if span <= eps
    y = x;
    return;
end
step = span / max(levels - 1, 1);
y = round((x - min(x)) / step) * step + min(x);
end

function value = rms_value(x)
value = sqrt(mean(x(:) .^ 2));
end

function textOut = safe_name(textIn)
textOut = regexprep(char(textIn), '[^A-Za-z0-9_]+', '_');
textOut = regexprep(textOut, '_+', '_');
textOut = regexprep(textOut, '^_|_$', '');
end
