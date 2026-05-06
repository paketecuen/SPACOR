function campaign = run_threewire_degradation(config, opts)
%RUN_THREEWIRE_DEGRADATION Canonical three-wire degradation campaign.

spacor.core.add_legacy_paths();
if nargin < 1 || isempty(config)
    config = spacor.threewire.default_config();
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
config.harmonicCandidates = spacor.core.get_option(config, ...
    'harmonicCandidates', [1 3 5 7]);
config.resultsDir = spacor.core.get_option(config, 'resultsDir', ...
    fullfile(spacor.core.lab_root(), 'results', 'threewire'));

outDir = spacor.core.get_option(opts, 'outputDir', ...
    fullfile(config.resultsDir, ['degradation_' mode]));
outDir = spacor.io.ensure_dir(outDir);
figureDir = spacor.io.ensure_dir(fullfile(outDir, 'figures'));
makeFigures = spacor.core.get_option(opts, 'makeFigures', true);

cases = threewire_cases(mode);
degradations = degradation_specs(mode);
nRows = numel(cases) * numel(degradations);
summaryRows = cell(nRows, 20);
parameterTables = cell(nRows, 1);
stabilityRows = cell(nRows, 13);
informationRows = cell(nRows, 14);
rankingTables = cell(nRows, 1);
candidateTables = cell(nRows, 1);
reports = cell(nRows, 1);
dataCases = cell(nRows, 1);

row = 0;
for c = 1:numel(cases)
    baseData = make_threewire_case(cases(c), config);
    for d = 1:numel(degradations)
        row = row + 1;
        deg = degradations(d);
        rng(9300 + 100 * c + d);
        data = apply_threewire_degradation(baseData, deg);
        data.label = sprintf('%s | %s', cases(c).name, deg.name);
        report = spacor.threewire.characterize(data, config, ...
            struct('candidateMode', 'mixed_blackbox', ...
            'windowCycles', config.windowCycles, ...
            'hopCycles', config.hopCycles));
        reports{row} = report;
        dataCases{row} = data;

        pass = is_threewire_degradation_pass(report, deg);
        summaryRows(row, :) = {row, cases(c).name, deg.name, ...
            cases(c).truthModel, report.canonical.status, ...
            report.canonical.selectedModel, report.recommendation.topology, ...
            report.recommendation.terms, report.canonical.family, pass, ...
            report.recommendation.selectedCoverage, ...
            report.recommendation.adequateCoverage, ...
            report.recommendation.dominantFraction, ...
            report.recommendation.medianCondition, ...
            report.recommendation.medianPowerResidual, ...
            report.recommendation.medianEnergyResidual, ...
            report.input.voltageKvlResidual, report.input.currentKclResidual, ...
            deg.noiseSnrDb, deg.operation};

        parameterTables{row} = threewire_parameter_table(row, cases(c), ...
            deg, report);
        stabilityRows(row, :) = {row, cases(c).name, deg.name, pass, ...
            report.recommendation.selectedCoverage, ...
            report.recommendation.adequateCoverage, ...
            report.recommendation.medianCondition, ...
            report.recommendation.medianSigmaMin, ...
            report.recommendation.medianPowerResidual, ...
            report.recommendation.medianEnergyResidual, ...
            report.input.voltageKvlResidual, report.input.currentKclResidual, ...
            deg.noiseSnrDb};
        informationRows(row, :) = {row, cases(c).name, deg.name, ...
            data.fs, data.f0, numel(data.t), data.t(end) - data.t(1), ...
            rms_value(data.vab), rms_value(data.vbc), rms_value(data.vca), ...
            rms_value(data.ia), rms_value(data.ib), rms_value(data.ic), ...
            deg.operation};

        ranking = report.ranking;
        ranking.CaseIndex = repmat(row, height(ranking), 1);
        ranking.CaseName = repmat({cases(c).name}, height(ranking), 1);
        ranking.Degradation = repmat({deg.name}, height(ranking), 1);
        ranking = movevars(ranking, {'CaseIndex', 'CaseName', ...
            'Degradation'}, 'Before', 1);
        rankingTables{row} = ranking;

        candidates = report.candidates;
        candidates.CaseIndex = repmat(row, height(candidates), 1);
        candidates.CaseName = repmat({cases(c).name}, height(candidates), 1);
        candidates.Degradation = repmat({deg.name}, height(candidates), 1);
        candidates = movevars(candidates, {'CaseIndex', 'CaseName', ...
            'Degradation'}, 'Before', 1);
        candidateTables{row} = candidates;

        if makeFigures
            safe = sprintf('%03d_%s_%s_parameters.png', row, ...
                safe_name(cases(c).name), safe_name(deg.name));
            try
                spacor.plot.parameter_evolution(report, data, ...
                    fullfile(figureDir, safe), ...
                    struct('truthMargin', config.plotTruthMargin));
            catch err
                warning('spacor:plotFailed', ...
                    'Three-wire degradation plot failed for %s/%s: %s', ...
                    cases(c).name, deg.name, err.message);
            end
        end
    end
end

summary = cell2table(summaryRows, 'VariableNames', { ...
    'CaseIndex', 'CaseName', 'Degradation', 'TruthModel', 'Status', ...
    'SelectedModel', 'SelectedTopology', 'SelectedTerms', ...
    'SelectedFamily', 'Pass', 'SelectedCoverage', 'AdequateCoverage', ...
    'DominantFraction', 'MedianCondition', 'MedianPowerResidual', ...
    'MedianEnergyResidual', 'VoltageKvlResidual', 'CurrentKclResidual', ...
    'NoiseSnrDb', 'Operation'});
parameters = vertcat(parameterTables{:});
stability = cell2table(stabilityRows, 'VariableNames', { ...
    'CaseIndex', 'CaseName', 'Degradation', 'Pass', ...
    'SelectedCoverage', 'AdequateCoverage', 'MedianCondition', ...
    'MedianSigmaMin', 'MedianPowerResidual', 'MedianEnergyResidual', ...
    'VoltageKvlResidual', 'CurrentKclResidual', 'NoiseSnrDb'});
information = cell2table(informationRows, 'VariableNames', { ...
    'CaseIndex', 'CaseName', 'Degradation', 'Fs', 'F0', 'Samples', ...
    'Duration', 'VabRms', 'VbcRms', 'VcaRms', 'IaRms', 'IbRms', ...
    'IcRms', 'Operation'});
method_ranking = vertcat(rankingTables{:});
candidates = vertcat(candidateTables{:});
summaryFigure = fullfile(figureDir, 'campaign_summary.png');
spacor.plot.campaign_summary(struct('summary', summary, ...
    'stability', stability, 'outputDir', outDir), summaryFigure);

writetable(summary, fullfile(outDir, 'summary.csv'));
writetable(parameters, fullfile(outDir, 'parameters.csv'));
writetable(stability, fullfile(outDir, 'stability.csv'));
writetable(information, fullfile(outDir, 'information.csv'));
writetable(method_ranking, fullfile(outDir, 'method_ranking.csv'));
writetable(candidates, fullfile(outDir, 'candidates.csv'));
save(fullfile(outDir, 'report.mat'), 'summary', 'parameters', ...
    'stability', 'information', 'method_ranking', 'candidates', ...
    'reports', 'dataCases', 'config', 'opts', 'mode', 'summaryFigure');

campaign = struct('summary', summary, 'parameters', parameters, ...
    'stability', stability, 'information', information, ...
    'method_ranking', method_ranking, 'candidates', candidates, ...
    'reports', {reports}, 'dataCases', {dataCases}, 'config', config, ...
    'mode', mode, 'outputDir', outDir, 'figureDir', figureDir, ...
    'summaryFigure', summaryFigure, ...
    'api', struct('name', 'spacor.campaigns.run_threewire_degradation'));
end

function cases = threewire_cases(mode)
cases = struct('name', {}, 'truthModel', {}, 'kind', {});
cases(end + 1) = struct('name', 'delta_g_unbalanced_multisine', ...
    'truthModel', 'delta_parallel_g', 'kind', 'delta_g');
cases(end + 1) = struct('name', 'wye_r_unbalanced_multisine', ...
    'truthModel', 'wye_series_r', 'kind', 'wye_r');
if strcmp(mode, 'smoke')
    cases = cases(1);
end
end

function specs = degradation_specs(mode)
specs = struct('name', {}, 'noiseSnrDb', {}, 'operation', {});
specs(end + 1) = make_deg('clean', Inf, 'none');
specs(end + 1) = make_deg('noise_60db', 60, 'noise');
if ~strcmp(mode, 'smoke')
    specs(end + 1) = make_deg('noise_40db', 40, 'noise');
    specs(end + 1) = make_deg('gain_mismatch_1pct', Inf, 'gain_mismatch');
    specs(end + 1) = make_deg('offset_1pct', Inf, 'offset');
    specs(end + 1) = make_deg('quantized_12bit', Inf, 'quantized');
    specs(end + 1) = make_deg('moving_average_5', Inf, 'moving_average');
end
end

function deg = make_deg(name, noiseSnrDb, operation)
deg = struct('name', name, 'noiseSnrDb', noiseSnrDb, ...
    'operation', operation);
end

function data = make_threewire_case(spec, config)
t = (0:round(config.duration * config.fs) - 1).' / config.fs;
w0 = 2 * pi * config.f0;
[va, vb, vc] = three_phase_voltage(t, w0);
vab = va - vb;
vbc = vb - vc;
vca = vc - va;

switch spec.kind
    case 'delta_g'
        truth = struct('Gab', 0.08, 'Gbc', 0.10, 'Gca', 0.12);
        iab = truth.Gab * vab;
        ibc = truth.Gbc * vbc;
        ica = truth.Gca * vca;
        ia = iab - ica;
        ib = ibc - iab;
        ic = ica - ibc;
    case 'wye_r'
        truth = struct('Ra', 12, 'Rb', 16, 'Rc', 20);
        ya = 1 / truth.Ra;
        yb = 1 / truth.Rb;
        yc = 1 / truth.Rc;
        vn = (ya * va + yb * vb + yc * vc) / (ya + yb + yc);
        ia = ya * (va - vn);
        ib = yb * (vb - vn);
        ic = yc * (vc - vn);
    otherwise
        error('Unknown three-wire degradation case "%s".', spec.kind);
end

data = struct('t', t, 'vab', vab, 'vbc', vbc, 'vca', vca, ...
    'ia', ia, 'ib', ib, 'ic', ic, 'fs', config.fs, ...
    'f0', config.f0, 'label', spec.name, 'truth', truth, ...
    'truthModel', spec.truthModel);
end

function data = apply_threewire_degradation(data, deg)
if isfinite(deg.noiseSnrDb)
    fields = {'vab', 'vbc', 'vca', 'ia', 'ib', 'ic'};
    for idx = 1:numel(fields)
        data.(fields{idx}) = add_noise(data.(fields{idx}), deg.noiseSnrDb);
    end
end
switch deg.operation
    case 'gain_mismatch'
        gains = struct('vab', 1.010, 'vbc', 0.990, 'vca', 1.004, ...
            'ia', 0.992, 'ib', 1.008, 'ic', 1.000);
        data = scale_fields(data, gains);
    case 'offset'
        fields = {'vab', 'vbc', 'vca', 'ia', 'ib', 'ic'};
        for idx = 1:numel(fields)
            f = fields{idx};
            data.(f) = data.(f) + 0.01 * rms_value(data.(f));
        end
    case 'quantized'
        fields = {'vab', 'vbc', 'vca', 'ia', 'ib', 'ic'};
        for idx = 1:numel(fields)
            data.(fields{idx}) = quantize_signal(data.(fields{idx}), 12);
        end
    case 'moving_average'
        fields = {'vab', 'vbc', 'vca', 'ia', 'ib', 'ic'};
        for idx = 1:numel(fields)
            data.(fields{idx}) = movmean(data.(fields{idx}), 5, ...
                'Endpoints', 'shrink');
        end
end
data.noiseSnrDb = deg.noiseSnrDb;
data.degradation = deg;
data.power = struct('measured', data.vab(:) .* data.ia(:) + ...
    data.vbc(:) .* (data.ia(:) + data.ib(:)));
end

function data = scale_fields(data, gains)
fields = fieldnames(gains);
for idx = 1:numel(fields)
    data.(fields{idx}) = gains.(fields{idx}) * data.(fields{idx});
end
end

function pass = is_threewire_degradation_pass(report, deg)
status = string(report.canonical.status);
if any(strcmp(deg.operation, {'gain_mismatch', 'offset'}))
    pass = any(status == ["accepted", "measurement_inconsistent", ...
        "low_coverage"]);
else
    pass = any(status == ["accepted", "ambiguous_equivalent_family"]) && ...
        report.recommendation.selectedCoverage >= 0.5;
end
end

function tableOut = threewire_parameter_table(caseIndex, spec, deg, report)
selected = report.windows.Selected > 0;
if any(selected)
    model = string(report.windows.SelectedModel(selected));
    topology = string(report.windows.SelectedTopology(selected));
    terms = string(report.windows.SelectedTerms(selected));
    nParams = report.windows.SelectedN(selected);
    condition = report.windows.ConditionNumber(selected);
    powerResidual = report.windows.PowerResidual(selected);
else
    model = strings(0, 1);
    topology = strings(0, 1);
    terms = strings(0, 1);
    nParams = [];
    condition = [];
    powerResidual = [];
end
tableOut = table(caseIndex, string(spec.name), string(deg.name), ...
    string(spec.truthModel), string(report.canonical.selectedModel), ...
    mode_string(model), mode_string(topology), mode_string(terms), ...
    median_finite(nParams), median_finite(condition), ...
    median_finite(powerResidual), ...
    'VariableNames', {'CaseIndex', 'CaseName', 'Degradation', ...
    'TruthModel', 'SelectedModel', 'MedianWindowModel', ...
    'MedianTopology', 'MedianTerms', 'MedianNParameters', ...
    'MedianCondition', 'MedianPowerResidual'});
end

function [va, vb, vc] = three_phase_voltage(t, w0)
va = phase_multisine(t, w0, 0, [120 7 4], [1 3 5], [0 0.2 -0.3]);
vb = phase_multisine(t, w0, -2 * pi / 3, [118 5 3], ...
    [1 3 7], [0.05 -0.4 0.6]);
vc = phase_multisine(t, w0, 2 * pi / 3, [122 6 5], ...
    [1 5 7], [-0.04 0.3 -0.5]);
end

function y = phase_multisine(t, w0, shift, amps, harmonics, phases)
y = zeros(size(t));
for idx = 1:numel(harmonics)
    omega = harmonics(idx) * w0;
    y = y + amps(idx) * sin(omega * t + harmonics(idx) * shift + phases(idx));
end
end

function y = add_noise(x, snrDb)
sigma = rms_value(x) * 10 ^ (-snrDb / 20);
y = x(:) + sigma * randn(size(x(:)));
end

function y = quantize_signal(x, bits)
x = x(:);
levels = 2 ^ bits;
span = max(x) - min(x);
if span <= eps
    y = x;
    return;
end
step = span / max(levels - 1, 1);
y = round((x - min(x)) / step) * step + min(x);
end

function value = mode_string(values)
if isempty(values)
    value = "";
    return;
end
[u, ~, idx] = unique(values, 'stable');
counts = accumarray(idx, 1);
[~, maxIdx] = max(counts);
value = u(maxIdx);
end

function value = median_finite(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = median(x);
end
end

function value = rms_value(x)
value = sqrt(mean(x(:) .^ 2));
end

function textOut = safe_name(textIn)
textOut = regexprep(char(textIn), '[^A-Za-z0-9_]+', '_');
textOut = regexprep(textOut, '_+', '_');
textOut = regexprep(textOut, '^_|_$', '');
end
