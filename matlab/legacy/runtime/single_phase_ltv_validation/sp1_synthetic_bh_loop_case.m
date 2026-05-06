function data = sp1_synthetic_bh_loop_case(config, opts)
%SP1_SYNTHETIC_BH_LOOP_CASE Stationary flux-current B-H loop.
%
% The magnetic state is represented by a periodic flux linkage lambda(t).
% The current is generated from an anhysteretic nonlinear curve plus an
% optional branch-dependent hysteresis term:
%
%   i = i_anh(lambda) + i_hyst(lambda, sign(dlambda/dt)).
%
% The terminal voltage is:
%
%   v = R i + d lambda / dt.
%
% The area of the lambda-i loop is the magnetic energy lost per cycle.

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
fluxKind = char(nv3_get_option(opts, 'fluxKind', 'sine'));
lambdaComponents = nv3_get_option(opts, 'lambdaComponents', ...
    default_flux_components(fluxKind, truth));

lambdaEval = nv3_harmonic_eval(lambdaComponents, t, f0, 1);
lambda = lambdaEval.values{1};
dlambda = lambdaEval.values{2};
[i, iAnh, iHyst] = bh_current(lambda, dlambda, truth);
v = truth.R * i + dlambda;

noiseSnrDb = nv3_get_option(opts, 'noiseSnrDb', config.noiseSnrDb);
seed = nv3_get_option(opts, 'randomSeed', config.randomSeed);
[vN, iN] = add_noise(seed, noiseSnrDb, v, i);

data = struct();
data.model = ['bh_loop_' fluxKind];
data.label = ['Single-phase B-H loop | ' fluxKind];
data.fs = fs;
data.f0 = f0;
data.t = t;
data.v = vN;
data.i = iN;
data.clean = struct('v', v, 'i', i, 'iAnhysteretic', iAnh, ...
    'iHysteresis', iHyst, 'lambda', lambda, 'dlambda', dlambda);
data.truth = truth;
data.truth.law = 'bh_loop';
data.noiseSnrDb = noiseSnrDb;
data.fluxKind = fluxKind;
data.lambdaComponents = lambdaComponents;
data.harmonics = unique(lambdaComponents(:, 2)).';
data.power = struct('measured', vN(:) .* iN(:), 'clean', v(:) .* i(:), ...
    'copper', truth.R * i(:) .^ 2, 'magnetic', dlambda(:) .* i(:));
data.energy = energy_metrics(t, lambda, dlambda, i, truth.R, f0);

end

function truth = default_truth()
truth = struct();
truth.R = 0.20;
truth.lambdaPeak = 0.80;
truth.L0 = 0.080;
truth.saturationCubic = 12.0;
truth.hysteresisCurrent = 0.90;
truth.hysteresisSmoothing = 0.04;
end

function components = default_flux_components(kind, truth)
switch kind
    case 'sine'
        components = [truth.lambdaPeak, 1, -pi / 2];
    case 'distorted'
        components = [ ...
            truth.lambdaPeak,        1, -pi / 2; ...
            0.06 * truth.lambdaPeak, 3,  0.30; ...
            0.03 * truth.lambdaPeak, 5, -1.10];
    otherwise
        error('Unknown flux kind "%s".', kind);
end
end

function [i, iAnh, iHyst] = bh_current(lambda, dlambda, truth)
lambdaScale = max(abs(lambda));
lambdaScale = max(lambdaScale, eps);
u = lambda ./ lambdaScale;
iAnh = lambda ./ truth.L0 + truth.saturationCubic .* lambda .^ 3;
branch = tanh(dlambda ./ ...
    max(truth.hysteresisSmoothing * max(abs(dlambda)), eps));
shape = 0.25 + 0.75 * abs(u) .^ 0.5;
iHyst = truth.hysteresisCurrent .* shape .* branch;
i = iAnh + iHyst;
end

function metrics = energy_metrics(t, lambda, dlambda, i, R, f0)
period = median(diff(find_zero_crossings(t, lambda)));
if ~isfinite(period) || period <= 0
    period = 1 / f0;
end
cycleSamples = max(4, round(period / median(diff(t))));
nCycles = floor(numel(t) / cycleSamples);
magneticEnergy = zeros(nCycles, 1);
copperEnergy = zeros(nCycles, 1);
for idx = 1:nCycles
    range = (idx - 1) * cycleSamples + (1:cycleSamples);
    magneticEnergy(idx) = trapz(t(range), i(range) .* dlambda(range));
    copperEnergy(idx) = trapz(t(range), R * i(range) .^ 2);
end
metrics = struct('MagneticLossPerCycle', magneticEnergy, ...
    'CopperLossPerCycle', copperEnergy, ...
    'MedianMagneticLossPerCycle', median(magneticEnergy, 'omitnan'), ...
    'MedianCopperLossPerCycle', median(copperEnergy, 'omitnan'));
end

function crossings = find_zero_crossings(t, x)
mask = x(1:end - 1) <= 0 & x(2:end) > 0;
idx = (1:numel(x) - 1).';
idx = idx(mask);
crossings = t(idx);
if numel(crossings) < 2
    crossings = [t(1); t(end)];
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
