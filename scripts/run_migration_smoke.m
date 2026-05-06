function results = run_migration_smoke()
%RUN_MIGRATION_SMOKE Verify canonical package facades on small synthetic cases.

root = startup_spacor_lab();
results = struct();

results.modelRegistry = run_model_registry_smoke();
results.derivatives = run_derivative_smoke();
results.singlephase = run_singlephase_smoke();
results.threewire = run_threewire_smoke();
results.singlephaseCampaign = run_singlephase_campaign_smoke();
results.threewireCampaign = run_threewire_campaign_smoke();

fprintf('Public API smoke passed.\n');
fprintf('  registry:    %d models, query=%s\n', ...
    results.modelRegistry.nModels, results.modelRegistry.queryModel);
fprintf('  derivative:  method=%s, gain=%.6g\n', ...
    results.derivatives.method, results.derivatives.gain);
fprintf('  singlephase: %s, selected=%s\n', ...
    results.singlephase.status, results.singlephase.selected);
fprintf('  threewire:   %s, selected=%s\n', ...
    results.threewire.status, results.threewire.selected);
fprintf('  sp1 campaign:%d cases -> %s\n', ...
    results.singlephaseCampaign.nCases, results.singlephaseCampaign.outputDir);
fprintf('  3w campaign: %d cases -> %s\n', ...
    results.threewireCampaign.nCases, results.threewireCampaign.outputDir);
fprintf('  root:        %s\n', root);
end

function summary = run_model_registry_smoke()
models = spacor.core.model_registry();
rlc = spacor.core.model_registry('r||l||c');

assert(numel(models) >= 8);
assert(strcmp(rlc.model, 'parallel_rlc'));
assert(isfield(rlc, 'minimumInformation'));

summary = struct();
summary.nModels = numel(models);
summary.queryModel = rlc.model;
end

function summary = run_derivative_smoke()
fs = 10000;
f0 = 50;
t = (0:fs / f0 * 4 - 1).' / fs;
w0 = 2 * pi * f0;
x = sin(w0 * t) + 0.05 * sin(3 * w0 * t + 0.2);
dx = w0 * cos(w0 * t) + 0.05 * 3 * w0 * cos(3 * w0 * t + 0.2);
opts = struct('f0', f0, 'harmonics', [1 3], ...
    'reference', {{x, dx}});
result = spacor.signal.estimate_derivatives(x, fs, [0 1], ...
    'harmonic', opts);

assert(isfield(result, 'api'));
assert(strcmp(result.api.name, 'spacor.signal.estimate_derivatives'));
assert(result.diagnostics.NRMSE(result.diagnostics.Order == 1) < 1e-10);

summary = struct();
summary.method = result.method;
summary.gain = result.diagnostics.Gain(result.diagnostics.Order == 1);
end

function summary = run_singlephase_smoke()
config = spacor.singlephase.default_config();
config.fs = 10000;
config.f0 = 50;
config.duration = 0.20;
config.windowCycles = 0.5;
config.hopCycles = 0.10;
config.harmonics = [1 3 5 7];

t = (0:round(config.duration * config.fs) - 1).' / config.fs;
w0 = 2 * pi * config.f0;
R = 10;
L = 50e-3;
harmonics = [1 3 5 7];
amps = [120 9 5 3];
phases = [0 0.4 -0.7 1.1];
v = zeros(size(t));
i = zeros(size(t));
for idx = 1:numel(harmonics)
    h = harmonics(idx);
    omega = h * w0;
    z = R + 1i * omega * L;
    v = v + amps(idx) * sin(omega * t + phases(idx));
    i = i + amps(idx) / abs(z) * sin(omega * t + phases(idx) - angle(z));
end

data = struct('t', t, 'v', v, 'i', i, 'fs', config.fs, ...
    'f0', config.f0, 'label', 'migration smoke single-phase series RL');
opts = struct('candidateMode', 'lti');
report = spacor.singlephase.characterize(data, config, opts);

assert(isfield(report, 'api'));
assert(strcmp(report.api.name, 'spacor.singlephase.characterize'));
assert(height(report.candidates) > 0);

summary = struct();
summary.status = char(report.canonical.status);
summary.selected = char(report.canonical.selectedModel);
summary.report = report;
end

function summary = run_threewire_smoke()
config = spacor.threewire.default_config();
config.fs = 10000;
config.f0 = 50;
config.duration = 0.20;
config.windowCycles = 0.5;
config.hopCycles = 0.10;
config.harmonics = [1 3 5 7];
config.harmonicCandidates = [1 3 5 7];

t = (0:round(config.duration * config.fs) - 1).' / config.fs;
w0 = 2 * pi * config.f0;
vpk = 120;
va = vpk * sin(w0 * t) + 7 * sin(3 * w0 * t + 0.2);
vb = vpk * sin(w0 * t - 2 * pi / 3) + 5 * sin(5 * w0 * t - 0.3);
vc = vpk * sin(w0 * t + 2 * pi / 3) + 4 * sin(7 * w0 * t + 0.5);
vab = va - vb;
vbc = vb - vc;
vca = vc - va;

G = [0.08 0.10 0.12];
iab = G(1) * vab;
ibc = G(2) * vbc;
ica = G(3) * vca;
ia = iab - ica;
ib = ibc - iab;
ic = ica - ibc;

data = struct('t', t, 'vab', vab, 'vbc', vbc, 'vca', vca, ...
    'ia', ia, 'ib', ib, 'ic', ic, 'fs', config.fs, 'f0', config.f0, ...
    'label', 'migration smoke three-wire delta G');
opts = struct('candidateMode', 'mixed_blackbox');
report = spacor.threewire.characterize(data, config, opts);

assert(isfield(report, 'api'));
assert(strcmp(report.api.name, 'spacor.threewire.characterize'));
assert(height(report.candidates) > 0);

summary = struct();
summary.status = char(report.canonical.status);
summary.selected = char(report.canonical.selectedModel);
summary.report = report;
end

function summary = run_singlephase_campaign_smoke()
campaign = spacor.campaigns.run_singlephase_quick();
required = {'summary.csv', 'parameters.csv', 'stability.csv', ...
    'information.csv', 'method_ranking.csv', 'report.mat'};
assert_required_files(campaign.outputDir, required);
assert_png_count(campaign.figureDir, 3);
assert(isfile(fullfile(campaign.figureDir, 'campaign_summary.png')));
assert(height(campaign.summary) >= 3);
summary = struct('nCases', height(campaign.summary), ...
    'outputDir', campaign.outputDir);
end

function summary = run_threewire_campaign_smoke()
campaign = spacor.campaigns.run_threewire_quick();
required = {'summary.csv', 'parameters.csv', 'stability.csv', ...
    'information.csv', 'method_ranking.csv', 'candidates.csv', ...
    'report.mat'};
assert_required_files(campaign.outputDir, required);
assert_png_count(campaign.figureDir, 2);
assert(isfile(fullfile(campaign.figureDir, 'campaign_summary.png')));
assert(height(campaign.summary) >= 2);
summary = struct('nCases', height(campaign.summary), ...
    'outputDir', campaign.outputDir);
end

function assert_png_count(outDir, minCount)
files = dir(fullfile(outDir, '*.png'));
assert(numel(files) >= minCount, ...
    'Expected at least %d PNG files in %s.', minCount, outDir);
end

function assert_required_files(outDir, required)
for idx = 1:numel(required)
    assert(isfile(fullfile(outDir, required{idx})), ...
        'Missing campaign output file: %s', required{idx});
end
end
