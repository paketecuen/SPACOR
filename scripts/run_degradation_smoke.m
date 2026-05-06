function results = run_degradation_smoke()
%RUN_DEGRADATION_SMOKE Verify canonical degradation campaigns.

startup_spacor_lab();

results = struct();
results.singlephase = spacor.campaigns.run_singlephase_degradation([], ...
    struct('mode', 'smoke'));
results.threewire = spacor.campaigns.run_threewire_degradation([], ...
    struct('mode', 'smoke'));

assert_standard_outputs(results.singlephase.outputDir, false);
assert_standard_outputs(results.threewire.outputDir, true);
assert_png_count(results.singlephase.figureDir, 2);
assert_png_count(results.threewire.figureDir, 2);
assert(isfile(fullfile(results.singlephase.figureDir, 'campaign_summary.png')));
assert(isfile(fullfile(results.threewire.figureDir, 'campaign_summary.png')));

fprintf('Degradation smoke passed.\n');
fprintf('  singlephase: %d cases -> %s\n', ...
    height(results.singlephase.summary), results.singlephase.outputDir);
fprintf('  threewire:   %d cases -> %s\n', ...
    height(results.threewire.summary), results.threewire.outputDir);
end

function assert_standard_outputs(outDir, hasCandidates)
required = {'summary.csv', 'parameters.csv', 'stability.csv', ...
    'information.csv', 'method_ranking.csv', 'report.mat'};
if hasCandidates
    required{end + 1} = 'candidates.csv';
end
for idx = 1:numel(required)
    assert(isfile(fullfile(outDir, required{idx})), ...
        'Missing degradation output file: %s', required{idx});
end
end

function assert_png_count(outDir, minCount)
files = dir(fullfile(outDir, '*.png'));
assert(numel(files) >= minCount, ...
    'Expected at least %d PNG files in %s.', minCount, outDir);
end
