function data = sp1_synthetic_timevarying_l_case(config, opts)
%SP1_SYNTHETIC_TIMEVARYING_L_CASE Synthetic R + explicit L(t) inductor.
%
% User-paper model:
%   v = R i + d(L i)/dt = R i + L(t) di/dt + i dL/dt
%
% The default law is loss_like, where dL/dt is prescribed as a positive
% iron-loss-equivalent resistance Rfe(t).  A current_saturation law is also
% available to expose the sign-changing dL/dt challenge.

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
law = char(nv3_get_option(opts, 'law', 'loss_like'));
iComponents = nv3_get_option(opts, 'iComponents', ...
    [11.0, 1, -0.20; 1.40, 3, 0.65; 0.70, 5, -1.10]);

iEval = nv3_harmonic_eval(iComponents, t, f0, 1);
i = iEval.values{1};
di = iEval.values{2};
[L, dLdt, Lcoef, lambda] = ltv_law(law, i, di, t, truth);
v = truth.R * i + L .* di + dLdt .* i;

noiseSnrDb = nv3_get_option(opts, 'noiseSnrDb', config.noiseSnrDb);
seed = nv3_get_option(opts, 'randomSeed', config.randomSeed);
[vN, iN] = add_noise(seed, noiseSnrDb, v, i);

data = struct();
data.model = ['series_r_ltv_' law];
data.label = ['Single-phase series R + L(t) | ' law];
data.fs = fs;
data.f0 = f0;
data.t = t;
data.v = vN;
data.i = iN;
data.clean = struct('v', v, 'i', i, 'di', di, 'L', L, ...
    'Lcoef', Lcoef, 'Linc', Lcoef, 'dLdt', dLdt, ...
    'lambda', lambda, 'Rloss', dLdt, 'Reff', truth.R + dLdt);
data.truth = truth;
data.truth.law = law;
data.currentComponents = iComponents;
data.noiseSnrDb = noiseSnrDb;
data.power = struct('measured', vN(:) .* iN(:), 'clean', v(:) .* i(:));

end

function truth = default_truth()
truth = struct();
truth.R = 0.65;
truth.L0 = 0.080;
truth.RfeBase = 0.020;
truth.RfeGain = 0.060;
truth.Is = 7.0;
truth.Lmin = 0.035;
truth.Lsat0 = 0.180;
end

function [L, dLdt, Lcoef, lambda] = ltv_law(law, i, di, t, truth)
switch law
    case 'loss_like'
        dLdt = truth.RfeBase + truth.RfeGain * (i ./ truth.Is) .^ 2;
        L = truth.L0 + cumtrapz(t, dLdt);
        Lcoef = L;
        lambda = L .* i;
    case 'current_saturation'
        u = i ./ truth.Is;
        L = truth.Lmin + (truth.Lsat0 - truth.Lmin) ./ cosh(u) .^ 2;
        dLdi = -2 * (truth.Lsat0 - truth.Lmin) .* ...
            tanh(u) ./ (truth.Is .* cosh(u) .^ 2);
        dLdt = dLdi .* di;
        Lcoef = L + i .* dLdi;
        lambda = L .* i;
    otherwise
        error('Unknown LTV law "%s".', law);
end
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
