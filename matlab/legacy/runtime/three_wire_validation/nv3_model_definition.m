function def = nv3_model_definition(model)
%NV3_MODEL_DEFINITION Three-wire identifiable model definitions.

model = char(model);
aliases = containers.Map( ...
    {'delta_parallel_r', 'delta_parallel_l', 'delta_parallel_c', ...
    'delta_parallel_rl', 'delta_parallel_rc', 'delta_parallel_lc', ...
    'delta_parallel_rlc', 'wye_series_c'}, ...
    {'delta_parallel_g', 'delta_parallel_gamma', 'delta_parallel_c', ...
    'delta_parallel_gl', 'delta_parallel_gc', 'delta_parallel_gammac', ...
    'delta_parallel_glc', 'wye_series_gamma'});
if isKey(aliases, model)
    model = aliases(model);
end

def = struct();
def.model = model;

switch model
    case 'delta_parallel_g'
        def.topology = 'delta_parallel';
        def.terms = {'G'};
    case 'delta_parallel_gamma'
        def.topology = 'delta_parallel';
        def.terms = {'Gamma'};
    case 'delta_parallel_c'
        def.topology = 'delta_parallel';
        def.terms = {'C'};
    case 'delta_parallel_gl'
        def.topology = 'delta_parallel';
        def.terms = {'G', 'Gamma'};
    case 'delta_parallel_gc'
        def.topology = 'delta_parallel';
        def.terms = {'G', 'C'};
    case 'delta_parallel_gammac'
        def.topology = 'delta_parallel';
        def.terms = {'Gamma', 'C'};
    case 'delta_parallel_glc'
        def.topology = 'delta_parallel';
        def.terms = {'G', 'Gamma', 'C'};
    case 'wye_series_r'
        def.topology = 'wye_series';
        def.terms = {'R'};
    case 'wye_series_l'
        def.topology = 'wye_series';
        def.terms = {'L'};
    case 'wye_series_gamma'
        def.topology = 'wye_series';
        def.terms = {'Gamma'};
    case 'wye_series_rl'
        def.topology = 'wye_series';
        def.terms = {'R', 'L'};
    case 'wye_series_rc'
        def.topology = 'wye_series';
        def.terms = {'R', 'Gamma'};
    case 'wye_series_lc'
        def.topology = 'wye_series';
        def.terms = {'L', 'Gamma'};
    case 'wye_series_rlc'
        def.topology = 'wye_series';
        def.terms = {'R', 'L', 'Gamma'};
    otherwise
        error('Unknown three-wire model "%s".', model);
end

switch def.topology
    case 'delta_parallel'
        def.branches = {'ab', 'bc', 'ca'};
    case 'wye_series'
        def.branches = {'a', 'b', 'c'};
end
def.thetaNames = theta_names(def);
def.parameterNames = parameter_names(def);

end

function names = theta_names(def)
names = {};
for termIdx = 1:numel(def.terms)
    term = def.terms{termIdx};
    for branchIdx = 1:numel(def.branches)
        names{end + 1} = [term def.branches{branchIdx}]; %#ok<AGROW>
    end
end
end

function names = parameter_names(def)
names = def.thetaNames;
if strcmp(def.topology, 'delta_parallel')
    if any(strcmp(def.terms, 'G'))
        names = [names, {'Rab', 'Rbc', 'Rca'}];
    end
    if any(strcmp(def.terms, 'Gamma'))
        names = [names, {'Lab', 'Lbc', 'Lca'}];
    end
elseif strcmp(def.topology, 'wye_series')
    if any(strcmp(def.terms, 'Gamma'))
        names = [names, {'Ca', 'Cb', 'Cc'}];
    end
end
end
