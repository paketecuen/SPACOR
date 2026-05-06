function results = run_fixture_regression_tests()
%RUN_FIXTURE_REGRESSION_TESTS Regression tests using local synthetic fixtures.

root = startup_spacor_lab();
fixtureDir = fullfile(root, 'fixtures', 'synthetic');
manifest = spacor.fixtures.write_synthetic_fixtures(fixtureDir);

singlephase = load(fullfile(fixtureDir, 'singlephase_series_rl.mat'));
spReport = spacor.singlephase.characterize(singlephase.data, ...
    singlephase.config, struct('candidateMode', 'lti'));
assert(strcmp(spReport.api.name, 'spacor.singlephase.characterize'));
assert(strcmp(spReport.canonical.selectedModel, 'lti_series_rl'));
assert(startsWith(string(spReport.canonical.status), "accepted"));

threewire = load(fullfile(fixtureDir, 'threewire_delta_g.mat'));
twReport = spacor.threewire.characterize(threewire.data, ...
    threewire.config, struct('candidateMode', 'mixed_blackbox'));
assert(strcmp(twReport.api.name, 'spacor.threewire.characterize'));
assert(height(twReport.candidates) > 0);
assert(strlength(string(twReport.canonical.selectedModel)) > 0);

results = struct();
results.manifest = manifest;
results.singlephaseSelected = string(spReport.canonical.selectedModel);
results.threewireSelected = string(twReport.canonical.selectedModel);
results.fixtureDir = fixtureDir;

outDir = spacor.io.ensure_dir(fullfile(root, 'results', 'regression'));
save(fullfile(outDir, 'latest_fixture_regression.mat'), 'results', ...
    'spReport', 'twReport');

fprintf('Fixture regression passed.\n');
fprintf('  singlephase: %s\n', results.singlephaseSelected);
fprintf('  threewire:   %s\n', results.threewireSelected);
fprintf('  fixtures:    %s\n', fixtureDir);
end
