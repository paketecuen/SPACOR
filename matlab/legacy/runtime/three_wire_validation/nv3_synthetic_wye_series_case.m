function data = nv3_synthetic_wye_series_case(config, opts)
%NV3_SYNTHETIC_WYE_SERIES_CASE Three-wire wye series R/L/Gamma case.
%
% Branch model:
%   v_xn = R_x i_x + L_x di_x/dt + Gamma_x integral(i_x)dt,
% with i_a+i_b+i_c=0.  Gamma is the series elastance 1/C.
% Measured line voltages:
%   v_ab = v_an - v_bn, v_bc = v_bn - v_cn, v_ca = v_cn - v_an.

if nargin < 1 || isempty(config)
    config = nv3_default_config();
end
if nargin < 2
    opts = struct();
end

fs = nv3_get_option(opts, 'fs', config.fs);
f0 = nv3_get_option(opts, 'f0', config.f0);
duration = nv3_get_option(opts, 'duration', config.duration);
t = (0:round(duration * fs) - 1).' / fs;

model = char(nv3_get_option(opts, 'model', 'wye_series_rl'));
def = nv3_model_definition(model);
if ~strcmp(def.topology, 'wye_series')
    error('Model "%s" is not a wye series model.', model);
end

truth = complete_truth(nv3_get_option(opts, 'truth', default_truth()));
iaComponents = nv3_get_option(opts, 'iaComponents', ...
    [10.0, 1, 0.10; 0.90, 5, -0.30; 0.45, 7, 1.20]);
ibComponents = nv3_get_option(opts, 'ibComponents', ...
    [8.0, 1, -2.05; 0.70, 5, 0.85; 0.25, 7, -1.40]);

iaEval = nv3_harmonic_eval(iaComponents, t, f0, 1);
ibEval = nv3_harmonic_eval(ibComponents, t, f0, 1);
ia = iaEval.values{1};
ib = ibEval.values{1};
ic = -ia - ib;
dia = iaEval.values{2};
dib = ibEval.values{2};
dic = -dia - dib;
qia = iaEval.primitive;
qib = ibEval.primitive;
qic = -qia - qib;

van = phase_voltage(def, truth, 'a', ia, dia, qia);
vbn = phase_voltage(def, truth, 'b', ib, dib, qib);
vcn = phase_voltage(def, truth, 'c', ic, dic, qic);
vab = van - vbn;
vbc = vbn - vcn;
vca = vcn - van;

noiseSnrDb = nv3_get_option(opts, 'noiseSnrDb', config.noiseSnrDb);
seed = nv3_get_option(opts, 'randomSeed', config.randomSeed);
[vabN, vbcN, vcaN, iaN, ibN, icN] = add_noise(seed, noiseSnrDb, vab, vbc, vca, ia, ib, ic);

data = struct();
data.model = model;
data.label = ['Three-wire ' strrep(model, '_', ' ')];
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
    'ia', ia, 'ib', ib, 'ic', ic, 'dia', dia, 'dib', dib, 'dic', dic, ...
    'qia', qia, 'qib', qib, 'qic', qic, ...
    'van', van, 'vbn', vbn, 'vcn', vcn);
data.truth = truth;
data.power = struct('measured', three_wire_power(vabN, vbcN, iaN, ibN), ...
    'clean', three_wire_power(vab, vbc, ia, ib));
data.currentComponents = struct('ia', iaComponents, 'ib', ibComponents);
data.noiseSnrDb = noiseSnrDb;

end

function truth = default_truth()
truth = struct();
truth.Ra = 0.50;
truth.Rb = 1.20;
truth.Rc = 0.80;
truth.La = 5.0e-3;
truth.Lb = 8.0e-3;
truth.Lc = 3.0e-3;
truth.Gammaa = 1 / 450e-6;
truth.Gammab = 1 / 320e-6;
truth.Gammac = 1 / 680e-6;
truth.Ca = 1 / truth.Gammaa;
truth.Cb = 1 / truth.Gammab;
truth.Cc = 1 / truth.Gammac;
end

function truth = complete_truth(truth)
if ~isfield(truth, 'Gammaa') && isfield(truth, 'Ca')
    truth.Gammaa = 1 / truth.Ca;
end
if ~isfield(truth, 'Gammab') && isfield(truth, 'Cb')
    truth.Gammab = 1 / truth.Cb;
end
if ~isfield(truth, 'Gammac') && isfield(truth, 'Cc')
    truth.Gammac = 1 / truth.Cc;
end
if ~isfield(truth, 'Ca') && isfield(truth, 'Gammaa')
    truth.Ca = 1 / truth.Gammaa;
end
if ~isfield(truth, 'Cb') && isfield(truth, 'Gammab')
    truth.Cb = 1 / truth.Gammab;
end
if ~isfield(truth, 'Cc') && isfield(truth, 'Gammac')
    truth.Cc = 1 / truth.Gammac;
end
end

function v = phase_voltage(def, truth, branch, i, di, qi)
v = zeros(size(i));
for termIdx = 1:numel(def.terms)
    term = def.terms{termIdx};
    name = [term branch];
    switch term
        case 'R'
            v = v + truth.(name) * i;
        case 'L'
            v = v + truth.(name) * di;
        case 'Gamma'
            v = v + truth.(name) * qi;
        otherwise
            error('Unknown wye series term "%s".', term);
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
