function config = sp1_default_config()
%SP1_DEFAULT_CONFIG Defaults for single-phase LTV validation.

root = fileparts(mfilename('fullpath'));
config = struct();
config.fs = 20000;
config.f0 = 50;
config.duration = 0.5;
config.harmonics = [1 3 5 7 9 11 13];
config.windowCycles = 0.5;
config.hopCycles = 0.05;
config.conditionLimit = 1e7;
config.noiseSnrDb = Inf;
config.randomSeed = 23;
config.resultsDir = fullfile(root, 'results');
config.plotTruthMargin = 0.10;

end
