function config = default_config()
%DEFAULT_CONFIG Canonical defaults for three-wire SPACOR characterization.

spacor.core.add_legacy_paths();
config = nv3_default_config();
config.resultsDir = fullfile(spacor.core.lab_root(), ...
    'results', 'threewire');
end

