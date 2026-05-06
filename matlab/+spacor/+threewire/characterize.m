function report = characterize(data, config, opts)
%CHARACTERIZE Characterize a three-wire black-box load.
%
% Canonical package entry point for the current three-wire framework.  The
% public facade delegates to the validated research runtime while the package
% namespace is kept stable for reproducibility.

spacor.core.add_legacy_paths();
if nargin < 2 || isempty(config)
    config = spacor.threewire.default_config();
end
if nargin < 3
    opts = struct();
end

report = nv3_characterize_three_wire_load(data, config, opts);
report.api = struct( ...
    'name', 'spacor.threewire.characterize', ...
    'legacyFunction', 'nv3_characterize_three_wire_load', ...
    'releaseStage', 'public_facade');
report.canonical = struct( ...
    'status', report.recommendation.status, ...
    'selectedModel', report.recommendation.model, ...
    'family', report.recommendation.family);
end
