function models = model_registry(query, opts)
%MODEL_REGISTRY Metadata for supported SPACOR circuit models.
%
% The registry keeps the physical vocabulary in one place.  It intentionally
% preserves the legacy field names used by early prototypes while adding a
% small canonical layer for newer campaigns.

if nargin < 1
    query = [];
end
if nargin < 2
    opts = struct(); %#ok<NASGU>
end

items = build_models();

if nargin < 1 || isempty(query)
    models = items;
    return;
end

queryText = lower(char(query));
if any(strcmp(queryText, {'singlephase', 'single_phase', 'monofasica', 'all'}))
    models = items;
    return;
end

name = canonical_model(query);
idx = find(strcmp({items.model}, name), 1);
if isempty(idx)
    error('spacor:unknownCircuitModel', ...
        'Unknown circuit model "%s".', char(query));
end
models = items(idx);
end

function items = build_models()
items = struct('model', {}, 'label', {}, 'domain', {}, 'topology', {}, ...
    'parameterNames', {}, 'positiveParameterNames', {}, 'truth', {}, ...
    'minHarmonics', {}, 'minimumInformation', {}, 'stress', {}, ...
    'fitInitials', {}, 'constraints', {}, 'api', {});

items(end + 1) = make_model('series_rl', 'R+L', 'series', ...
    {'R', 'L'}, struct('R', 5.5, 'L', 5e-3), 1, false, ...
    [5.5, 5e-3; 2.0, 2e-3; 10, 10e-3]);

items(end + 1) = make_model('parallel_rl', 'R||L', 'parallel', ...
    {'G', 'R', 'Gamma', 'L'}, ...
    struct('G', 1 / 5.5, 'R', 5.5, 'Gamma', 1 / 5e-3, 'L', 5e-3), ...
    1, false, [1 / 5.5, 5e-3; 0.5, 2e-3; 0.05, 10e-3]);

items(end + 1) = make_model('series_rc', 'R+C', 'series', ...
    {'R', 'C'}, struct('R', 5.5, 'C', 250e-6), 1, false, ...
    [5.5, 250e-6; 2.0, 100e-6; 10, 1e-3]);

items(end + 1) = make_model('parallel_rc', 'R||C', 'parallel', ...
    {'G', 'R', 'C'}, struct('G', 1 / 5.5, 'R', 5.5, 'C', 250e-6), ...
    1, false, [1 / 5.5, 250e-6; 0.5, 100e-6; 0.05, 1e-3]);

items(end + 1) = make_model('series_lc', 'L+C', 'series', ...
    {'L', 'C'}, struct('L', 5e-3, 'C', 50e-6), 2, false, ...
    [5e-3, 50e-6; 2e-3, 100e-6; 10e-3, 25e-6]);

items(end + 1) = make_model('parallel_lc', 'L||C', 'parallel', ...
    {'Gamma', 'L', 'C'}, struct('Gamma', 1 / 5e-3, 'L', 5e-3, ...
    'C', 50e-6), 2, false, ...
    [5e-3, 50e-6; 2e-3, 100e-6; 10e-3, 25e-6]);

items(end + 1) = make_model('r_plus_rl', 'R||(R+L)', ...
    'parallel_series', {'Gp', 'Rp', 'Rser', 'Lser'}, ...
    struct('Gp', 1 / 0.77, 'Rp', 0.77, 'Rser', 0.148, ...
    'Lser', 5e-3), 2, false, ...
    [1 / 0.77, 0.148, 5e-3; 1.0, 0.5, 1e-3; 0.05, 0.6, 1e-3]);

items(end + 1) = make_model('series_rlc', 'R+L+C stress', 'series', ...
    {'R', 'L', 'C', 'Gamma'}, ...
    struct('R', 5.5, 'L', 5e-3, 'C', 250e-6, ...
    'Gamma', 1 / 250e-6), 2, true, ...
    [5.5, 5e-3, 250e-6; 2.0, 2e-3, 100e-6; 10, 10e-3, 1e-3]);

items(end + 1) = make_model('parallel_rlc', 'R||L||C stress', ...
    'parallel', {'G', 'R', 'C', 'Gamma', 'L'}, ...
    struct('G', 1 / 5.5, 'R', 5.5, 'C', 250e-6, ...
    'Gamma', 1 / 5e-3, 'L', 5e-3), 2, true, ...
    [1 / 5.5, 250e-6, 1 / 5e-3; 0.5, 100e-6, 500; ...
    0.05, 1e-3, 50]);
end

function item = make_model(model, label, topology, parameterNames, truth, ...
    minHarmonics, stress, fitInitials)
item = struct();
item.model = model;
item.label = label;
item.domain = 'singlephase';
item.topology = topology;
item.parameterNames = parameterNames;
item.positiveParameterNames = parameterNames;
item.truth = truth;
item.minHarmonics = minHarmonics;
item.minimumInformation = struct( ...
    'minHarmonics', minHarmonics, ...
    'note', minimum_information_note(minHarmonics, stress));
item.stress = stress;
item.fitInitials = fitInitials;
item.constraints = struct('passive', true, ...
    'positiveParameterNames', {parameterNames});
item.api = struct('name', 'spacor.core.model_registry', ...
    'releaseStage', 'public_api');
end

function note = minimum_information_note(minHarmonics, stress)
if stress
    note = ['requires rich windows; pure fundamental is not a reliable ' ...
        'identification gate'];
elseif minHarmonics <= 1
    note = 'can be estimated from a single informative frequency window';
else
    note = 'requires more than one independent excitation component';
end
end

function model = canonical_model(model)
switch lower(char(model))
    case {'rl', 'series_rl', 'rl_series', 'r+l'}
        model = 'series_rl';
    case {'parallel_rl', 'rl_parallel', 'r||l'}
        model = 'parallel_rl';
    case {'series_rc', 'rc_series', 'r+c'}
        model = 'series_rc';
    case {'rc', 'parallel_rc', 'rc_parallel', 'r||c'}
        model = 'parallel_rc';
    case {'series_lc', 'lc_series', 'l+c'}
        model = 'series_lc';
    case {'parallel_lc', 'lc_parallel', 'l||c'}
        model = 'parallel_lc';
    case {'r+rl', 'r_plus_rl', 'parallel_r_plus_series_rl', 'r||(r+l)'}
        model = 'r_plus_rl';
    case {'series_rlc', 'rlc_series', 'r+l+c', 'series_r_l_c'}
        model = 'series_rlc';
    case {'rlc', 'parallel_rlc', 'rlc_parallel', 'r||l||c', ...
            'parallel_r_l_c'}
        model = 'parallel_rlc';
    otherwise
        model = lower(char(model));
end
end
