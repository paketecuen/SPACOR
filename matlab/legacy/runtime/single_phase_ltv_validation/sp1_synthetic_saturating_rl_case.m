function data = sp1_synthetic_saturating_rl_case(config, opts)
%SP1_SYNTHETIC_SATURATING_RL_CASE Synthetic R + saturating inductor.
%
% Constitutive model:
%   v = R i + d(lambda(i))/dt
%   lambda(i) = Lmin i + (L0-Lmin) Is tanh(i/Is)
%   Linc(i) = d(lambda)/di = Lmin + (L0-Lmin) sech(i/Is)^2

if nargin < 1 || isempty(config)
    config = sp1_default_config();
end
if nargin < 2
    opts = struct();
end

fs = nv3_get_option(opts, 'fs', config.fs);
f0 = nv3_get_option(opts, 'f0', config.f0);
duration = nv3_get_option(opts, 'duration', config.duration);
t = (0:round(duration * fs) - 1).' / fs;

truth = nv3_get_option(opts, 'truth', default_truth());
iComponents = nv3_get_option(opts, 'iComponents', ...
    [11.0, 1, -0.20; 1.40, 3, 0.65; 0.70, 5, -1.10]);

iEval = nv3_harmonic_eval(iComponents, t, f0, 1);
i = iEval.values{1};
di = iEval.values{2};
lambda = sp1_flux_curve(i, truth);
Linc = sp1_incremental_inductance(i, truth);
v = truth.R * i + Linc .* di;

noiseSnrDb = nv3_get_option(opts, 'noiseSnrDb', config.noiseSnrDb);
seed = nv3_get_option(opts, 'randomSeed', config.randomSeed);
[vN, iN] = add_noise(seed, noiseSnrDb, v, i);

data = struct();
data.model = 'series_r_lsat';
data.label = 'Single-phase series R + saturating L';
data.fs = fs;
data.f0 = f0;
data.t = t;
data.v = vN;
data.i = iN;
data.clean = struct('v', v, 'i', i, 'di', di, ...
    'lambda', lambda, 'Linc', Linc);
data.truth = truth;
data.currentComponents = iComponents;
data.noiseSnrDb = noiseSnrDb;
data.power = struct('measured', vN(:) .* iN(:), 'clean', v(:) .* i(:));

end

function truth = default_truth()
truth = struct();
truth.R = 0.65;
truth.L0 = 0.180;
truth.Lmin = 0.035;
truth.Is = 7.0;
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
