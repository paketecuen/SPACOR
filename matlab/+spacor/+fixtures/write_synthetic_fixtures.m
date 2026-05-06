function manifest = write_synthetic_fixtures(outDir)
%WRITE_SYNTHETIC_FIXTURES Write deterministic fixture MAT files.

if nargin < 1 || isempty(outDir)
    outDir = fullfile(spacor.core.lab_root(), 'fixtures', 'synthetic');
end
outDir = spacor.io.ensure_dir(outDir);

[data, config, truth] = spacor.fixtures.singlephase_series_rl();
singlephaseFile = fullfile(outDir, 'singlephase_series_rl.mat');
save(singlephaseFile, 'data', 'config', 'truth');

[data, config, truth] = spacor.fixtures.threewire_delta_g(); %#ok<ASGLU>
threewireFile = fullfile(outDir, 'threewire_delta_g.mat');
save(threewireFile, 'data', 'config', 'truth');

manifest = table( ...
    ["singlephase_series_rl"; "threewire_delta_g"], ...
    [string(singlephaseFile); string(threewireFile)], ...
    ["lti_series_rl"; "delta_resistive"], ...
    'VariableNames', {'Fixture', 'Path', 'TruthModel'});
writetable(manifest, fullfile(outDir, 'manifest.csv'));
end
