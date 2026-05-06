function results = run_regression_tests(opts)
%RUN_REGRESSION_TESTS Run the canonical SPACOR public regression smoke suite.
%
% This is intentionally a light regression suite. It verifies that the
% canonical API, synthetic fixtures, and synthetic degradation campaigns still
% execute and produce standard outputs.  This public suite intentionally does
% not require or access private measurement data.

if nargin < 1
    opts = struct();
end

root = startup_spacor_lab();
started = datetime('now');
timer = tic();

results = struct();
results.root = root;
results.started = started;
results.fixtures = run_fixture_regression_tests();
results.migration = run_migration_smoke();
results.degradation = run_degradation_smoke();
results.elapsedSeconds = toc(timer);
results.finished = datetime('now');

outDir = spacor.io.ensure_dir(fullfile(root, 'results', 'regression'));
save(fullfile(outDir, 'latest_regression.mat'), 'results');

if spacor.core.get_option(opts, 'writeTextSummary', true)
    write_text_summary(results, fullfile(outDir, 'latest_regression.txt'));
end

fprintf('SPACOR regression passed in %.1f s.\n', results.elapsedSeconds);
fprintf('  output: %s\n', outDir);
end

function write_text_summary(results, outFile)
fid = fopen(outFile, 'w');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'SPACOR Public Regression\n');
fprintf(fid, 'Started:  %s\n', string(results.started));
fprintf(fid, 'Finished: %s\n', string(results.finished));
fprintf(fid, 'Elapsed:  %.3f s\n\n', results.elapsedSeconds);
fprintf(fid, 'Fixture single-phase:   %s\n', ...
    results.fixtures.singlephaseSelected);
fprintf(fid, 'Fixture three-wire:     %s\n', ...
    results.fixtures.threewireSelected);
fprintf(fid, 'Migration single-phase: %s\n', ...
    results.migration.singlephase.selected);
fprintf(fid, 'Migration three-wire:   %s\n', ...
    results.migration.threewire.selected);
fprintf(fid, 'Degradation SP cases:   %d\n', ...
    height(results.degradation.singlephase.summary));
fprintf(fid, 'Degradation 3W cases:   %d\n', ...
    height(results.degradation.threewire.summary));
end
