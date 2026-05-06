function config = default_config()
%DEFAULT_CONFIG Canonical defaults for single-phase SPACOR characterization.

spacor.core.add_legacy_paths();
config = sp1_default_config();
config.resultsDir = fullfile(spacor.core.lab_root(), ...
    'results', 'singlephase');
end

