function config = nv3_default_config()
%NV3_DEFAULT_CONFIG Defaults for three-wire numerical validation.

config = struct();
config.fs = 20000;
config.f0 = 50;
config.duration = 0.6;

config.harmonics = [1 3 5 7];
config.harmonicCandidates = [1 3 5 7 9 11 13];
config.adaptiveHarmonics = true;
config.harmonicSelectionMode = 'energy';
config.activeThreshold = 1e-3;
config.maxBasisCondition = 1e10;

config.windowCycles = 0.5;
config.hopCycles = 0.1;
config.minRowsFactor = 4;
config.conditionLimit = 1e9;

config.noiseSnrDb = Inf;
config.randomSeed = 7;

config.plotTruthMargin = 0.10;
config.resultsDir = fullfile(fileparts(mfilename('fullpath')), 'results');

end
