function report = characterize(data, config, opts)
%CHARACTERIZE Characterize a single-phase black-box load.
%
% Canonical package entry point for the current single-phase framework.
% The public facade delegates to the validated research runtime while the
% package namespace is kept stable for reproducibility.

spacor.core.add_legacy_paths();
if nargin < 2 || isempty(config)
    config = spacor.singlephase.default_config();
end
if nargin < 3
    opts = struct();
end

report = sp1_characterize_physical_load(data, config, opts);
report.api = struct( ...
    'name', 'spacor.singlephase.characterize', ...
    'legacyFunction', 'sp1_characterize_physical_load', ...
    'releaseStage', 'public_facade');
report.canonical = struct( ...
    'status', report.recommendation.decision, ...
    'selectedModel', report.recommendation.model, ...
    'family', report.recommendation.family);
end
