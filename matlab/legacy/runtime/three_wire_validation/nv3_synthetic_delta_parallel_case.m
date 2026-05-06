function data = nv3_synthetic_delta_parallel_case(config, opts)
%NV3_SYNTHETIC_DELTA_PARALLEL_CASE Three-wire delta parallel G/Gamma case.
%
% Branch model:
%   i_xy = G_xy v_xy + Gamma_xy integral(v_xy) dt
% Line currents:
%   i_a = i_ab - i_ca, i_b = i_bc - i_ab, i_c = i_ca - i_bc.

if nargin < 1 || isempty(config)
    config = nv3_default_config();
end
if nargin < 2
    opts = struct();
end

model = char(nv3_get_option(opts, 'model', 'delta_parallel_gl'));
def = nv3_model_definition(model);
if ~strcmp(def.topology, 'delta_parallel')
    error('Model "%s" is not a delta parallel model.', model);
end

fs = nv3_get_option(opts, 'fs', config.fs);
f0 = nv3_get_option(opts, 'f0', config.f0);
duration = nv3_get_option(opts, 'duration', config.duration);
t = (0:round(duration * fs) - 1).' / fs;

truth = complete_truth(nv3_get_option(opts, 'truth', default_truth()));
vabComponents = nv3_get_option(opts, 'vabComponents', ...
    [100, 1, 0.00; 8, 5, 0.40; 4, 7, -1.00]);
vbcComponents = nv3_get_option(opts, 'vbcComponents', ...
    [92, 1, -2.18; 6, 5, -0.75; 3, 7, 1.15]);

vabEval = nv3_harmonic_eval(vabComponents, t, f0, 1);
vbcEval = nv3_harmonic_eval(vbcComponents, t, f0, 1);
vab = vabEval.values{1};
vbc = vbcEval.values{1};
vca = -vab - vbc;
dvab = vabEval.values{2};
dvbc = vbcEval.values{2};
dvca = -dvab - dvbc;
qab = vabEval.primitive;
qbc = vbcEval.primitive;
qca = -qab - qbc;

iab = branch_current(def, truth, 'ab', vab, qab, dvab);
ibc = branch_current(def, truth, 'bc', vbc, qbc, dvbc);
ica = branch_current(def, truth, 'ca', vca, qca, dvca);
ia = iab - ica;
ib = ibc - iab;
ic = ica - ibc;

noiseSnrDb = nv3_get_option(opts, 'noiseSnrDb', config.noiseSnrDb);
seed = nv3_get_option(opts, 'randomSeed', config.randomSeed);
[vabN, vbcN, vcaN, iaN, ibN, icN] = add_noise(seed, noiseSnrDb, vab, vbc, vca, ia, ib, ic);

data = struct();
data.model = def.model;
data.label = ['Three-wire ' strrep(def.model, '_', ' ')];
data.fs = fs;
data.f0 = f0;
data.t = t;
data.vab = vabN;
data.vbc = vbcN;
data.vca = vcaN;
data.ia = iaN;
data.ib = ibN;
data.ic = icN;
data.clean = struct('vab', vab, 'vbc', vbc, 'vca', vca, ...
    'ia', ia, 'ib', ib, 'ic', ic, 'iab', iab, 'ibc', ibc, 'ica', ica, ...
    'qab', qab, 'qbc', qbc, 'qca', qca, ...
    'dvab', dvab, 'dvbc', dvbc, 'dvca', dvca);
data.truth = truth;
data.power = struct('measured', three_wire_power(vabN, vbcN, iaN, ibN), ...
    'clean', three_wire_power(vab, vbc, ia, ib));
data.voltageComponents = struct('vab', vabComponents, 'vbc', vbcComponents);
data.noiseSnrDb = noiseSnrDb;

end

function truth = default_truth()
truth = struct();
truth.Gab = 1 / 20;
truth.Gbc = 1 / 30;
truth.Gca = 1 / 15;
truth.Lab = 0.080;
truth.Lbc = 0.120;
truth.Lca = 0.060;
truth.Gammaab = 1 / truth.Lab;
truth.Gammabc = 1 / truth.Lbc;
truth.Gammaca = 1 / truth.Lca;
truth.Rab = 1 / truth.Gab;
truth.Rbc = 1 / truth.Gbc;
truth.Rca = 1 / truth.Gca;
truth.Cab = 80e-6;
truth.Cbc = 120e-6;
truth.Cca = 60e-6;
end

function truth = complete_truth(truth)
if ~isfield(truth, 'Rab') && isfield(truth, 'Gab')
    truth.Rab = 1 / truth.Gab;
end
if ~isfield(truth, 'Gab') && isfield(truth, 'Rab')
    truth.Gab = 1 / truth.Rab;
end
if ~isfield(truth, 'Rbc') && isfield(truth, 'Gbc')
    truth.Rbc = 1 / truth.Gbc;
end
if ~isfield(truth, 'Gbc') && isfield(truth, 'Rbc')
    truth.Gbc = 1 / truth.Rbc;
end
if ~isfield(truth, 'Rca') && isfield(truth, 'Gca')
    truth.Rca = 1 / truth.Gca;
end
if ~isfield(truth, 'Gca') && isfield(truth, 'Rca')
    truth.Gca = 1 / truth.Rca;
end
if ~isfield(truth, 'Gammaab') && isfield(truth, 'Lab')
    truth.Gammaab = 1 / truth.Lab;
end
if ~isfield(truth, 'Gammabc') && isfield(truth, 'Lbc')
    truth.Gammabc = 1 / truth.Lbc;
end
if ~isfield(truth, 'Gammaca') && isfield(truth, 'Lca')
    truth.Gammaca = 1 / truth.Lca;
end
if ~isfield(truth, 'Lab') && isfield(truth, 'Gammaab')
    truth.Lab = 1 / truth.Gammaab;
end
if ~isfield(truth, 'Lbc') && isfield(truth, 'Gammabc')
    truth.Lbc = 1 / truth.Gammabc;
end
if ~isfield(truth, 'Lca') && isfield(truth, 'Gammaca')
    truth.Lca = 1 / truth.Gammaca;
end
end

function i = branch_current(def, truth, branch, v, q, dv)
i = zeros(size(v));
for termIdx = 1:numel(def.terms)
    term = def.terms{termIdx};
    name = [term branch];
    switch term
        case 'G'
            i = i + truth.(name) * v;
        case 'Gamma'
            i = i + truth.(name) * q;
        case 'C'
            i = i + truth.(name) * dv;
    end
end
end

function p = three_wire_power(vab, vbc, ia, ib)
p = vab(:) .* ia(:) + vbc(:) .* (ia(:) + ib(:));
end

function varargout = add_noise(seed, snrDb, varargin)
rng(seed);
varargout = cell(size(varargin));
for idx = 1:numel(varargin)
    x = varargin{idx};
    if isfinite(snrDb)
        sigma = sqrt(mean(x(:) .^ 2)) * 10 ^ (-snrDb / 20);
        varargout{idx} = x + sigma * randn(size(x));
    else
        varargout{idx} = x;
    end
end
end
