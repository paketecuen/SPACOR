function campaign = run_threewire_quick(config, opts)
%RUN_THREEWIRE_QUICK Canonical quick three-wire SPACOR campaign.

spacor.core.add_legacy_paths();
if nargin < 1 || isempty(config)
    config = spacor.threewire.default_config();
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
config.harmonicCandidates = spacor.core.get_option(config, ...
    'harmonicCandidates', [1 3 5 7]);
config.resultsDir = spacor.core.get_option(config, 'resultsDir', ...
    fullfile(spacor.core.lab_root(), 'results', 'threewire'));

outDir = spacor.core.get_option(opts, 'outputDir', ...
    fullfile(config.resultsDir, 'quick'));
outDir = spacor.io.ensure_dir(outDir);
figureDir = spacor.io.ensure_dir(fullfile(outDir, 'figures'));
candidateMode = spacor.core.get_option(opts, 'candidateMode', ...
    'mixed_blackbox');

caseDefs = threewire_cases();
summaryRows = cell(numel(caseDefs), 17);
parameterTables = cell(numel(caseDefs), 1);
stabilityRows = cell(numel(caseDefs), 10);
informationRows = cell(numel(caseDefs), 11);
rankingTables = cell(numel(caseDefs), 1);
candidateTables = cell(numel(caseDefs), 1);
reports = cell(numel(caseDefs), 1);
dataCases = cell(numel(caseDefs), 1);

for idx = 1:numel(caseDefs)
    spec = caseDefs(idx);
    data = make_threewire_case(spec, config);
    report = spacor.threewire.characterize(data, config, ...
        struct('candidateMode', candidateMode, ...
        'windowCycles', config.windowCycles, ...
        'hopCycles', config.hopCycles));
    reports{idx} = report;
    dataCases{idx} = data;
    spacor.plot.parameter_evolution(report, data, ...
        fullfile(figureDir, [safe_name(spec.name) '_parameters.png']), ...
        struct('truthMargin', config.plotTruthMargin));

    pass = is_threewire_pass(report);
    summaryRows(idx, :) = {idx, spec.name, spec.truthModel, ...
        report.canonical.status, report.canonical.selectedModel, ...
        report.recommendation.topology, report.recommendation.terms, ...
        report.canonical.family, pass, ...
        report.recommendation.selectedCoverage, ...
        report.recommendation.adequateCoverage, ...
        report.recommendation.dominantFraction, ...
        report.recommendation.medianCondition, ...
        report.recommendation.medianPowerResidual, ...
        report.recommendation.medianEnergyResidual, ...
        report.input.voltageKvlResidual, report.input.currentKclResidual};

    parameterTables{idx} = threewire_parameter_table(idx, spec, report);
    stabilityRows(idx, :) = {idx, spec.name, pass, ...
        report.recommendation.selectedCoverage, ...
        report.recommendation.adequateCoverage, ...
        report.recommendation.medianCondition, ...
        report.recommendation.medianSigmaMin, ...
        report.recommendation.medianPowerResidual, ...
        report.input.voltageKvlResidual, report.input.currentKclResidual};
    informationRows(idx, :) = {idx, spec.name, data.fs, data.f0, ...
        numel(data.t), data.t(end) - data.t(1), ...
        rms_value(data.vab), rms_value(data.vbc), rms_value(data.vca), ...
        rms_value(data.ia), rms_value(data.ib)};

    ranking = report.ranking;
    ranking.CaseIndex = repmat(idx, height(ranking), 1);
    ranking.CaseName = repmat({spec.name}, height(ranking), 1);
    ranking = movevars(ranking, {'CaseIndex', 'CaseName'}, 'Before', 1);
    rankingTables{idx} = ranking;

    candidates = report.candidates;
    candidates.CaseIndex = repmat(idx, height(candidates), 1);
    candidates.CaseName = repmat({spec.name}, height(candidates), 1);
    candidates = movevars(candidates, {'CaseIndex', 'CaseName'}, 'Before', 1);
    candidateTables{idx} = candidates;
end

summary = cell2table(summaryRows, 'VariableNames', { ...
    'CaseIndex', 'CaseName', 'TruthModel', 'Status', 'SelectedModel', ...
    'SelectedTopology', 'SelectedTerms', 'SelectedFamily', 'Pass', ...
    'SelectedCoverage', 'AdequateCoverage', 'DominantFraction', ...
    'MedianCondition', 'MedianPowerResidual', 'MedianEnergyResidual', ...
    'VoltageKvlResidual', 'CurrentKclResidual'});
parameters = vertcat(parameterTables{:});
stability = cell2table(stabilityRows, 'VariableNames', { ...
    'CaseIndex', 'CaseName', 'Pass', 'SelectedCoverage', ...
    'AdequateCoverage', 'MedianCondition', 'MedianSigmaMin', ...
    'MedianPowerResidual', 'VoltageKvlResidual', 'CurrentKclResidual'});
information = cell2table(informationRows, 'VariableNames', { ...
    'CaseIndex', 'CaseName', 'Fs', 'F0', 'Samples', 'Duration', ...
    'VabRms', 'VbcRms', 'VcaRms', 'IaRms', 'IbRms'});
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
    'reports', 'dataCases', 'config', 'opts', 'summaryFigure');

campaign = struct('summary', summary, 'parameters', parameters, ...
    'stability', stability, 'information', information, ...
    'method_ranking', method_ranking, 'candidates', candidates, ...
    'reports', {reports}, 'dataCases', {dataCases}, 'config', config, ...
    'outputDir', outDir, 'figureDir', figureDir, ...
    'summaryFigure', summaryFigure, ...
    'api', struct('name', 'spacor.campaigns.run_threewire_quick'));
end

function cases = threewire_cases()
cases = struct('name', {}, 'truthModel', {}, 'kind', {});
cases(end + 1) = struct('name', 'delta_g_unbalanced_multisine', ...
    'truthModel', 'delta_parallel_g', 'kind', 'delta_g');
cases(end + 1) = struct('name', 'wye_r_unbalanced_multisine', ...
    'truthModel', 'wye_series_r', 'kind', 'wye_r');
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
        error('Unknown three-wire quick case "%s".', spec.kind);
end

data = struct('t', t, 'vab', vab, 'vbc', vbc, 'vca', vca, ...
    'ia', ia, 'ib', ib, 'ic', ic, 'fs', config.fs, ...
    'f0', config.f0, 'label', spec.name, 'truth', truth, ...
    'truthModel', spec.truthModel);
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

function pass = is_threewire_pass(report)
pass = any(strcmp(report.canonical.status, {'accepted', ...
    'ambiguous_equivalent_family'})) && ...
    report.recommendation.selectedCoverage >= 0.5 && ...
    report.input.voltageKvlResidual < 1e-10 && ...
    report.input.currentKclResidual < 1e-10;
end

function tableOut = threewire_parameter_table(caseIndex, spec, report)
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
tableOut = table(caseIndex, string(spec.name), string(spec.truthModel), ...
    string(report.canonical.selectedModel), mode_string(model), ...
    mode_string(topology), mode_string(terms), median_finite(nParams), ...
    median_finite(condition), median_finite(powerResidual), ...
    'VariableNames', {'CaseIndex', 'CaseName', 'TruthModel', ...
    'SelectedModel', 'MedianWindowModel', 'MedianTopology', ...
    'MedianTerms', 'MedianNParameters', 'MedianCondition', ...
    'MedianPowerResidual'});
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
