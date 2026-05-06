function [data, config, truth] = threewire_delta_g(config)
%THREEWIRE_DELTA_G Deterministic three-wire delta conductance fixture.

if nargin < 1 || isempty(config)
    config = spacor.threewire.default_config();
end

config.fs = spacor.core.get_option(config, 'fs', 10000);
config.f0 = spacor.core.get_option(config, 'f0', 50);
config.duration = spacor.core.get_option(config, 'duration', 0.20);
config.windowCycles = spacor.core.get_option(config, 'windowCycles', 0.5);
config.hopCycles = spacor.core.get_option(config, 'hopCycles', 0.10);
config.harmonics = spacor.core.get_option(config, 'harmonics', [1 3 5 7]);
config.harmonicCandidates = spacor.core.get_option(config, ...
    'harmonicCandidates', [1 3 5 7]);

t = (0:round(config.duration * config.fs) - 1).' / config.fs;
w0 = 2 * pi * config.f0;
vpk = 120;
va = vpk * sin(w0 * t) + 7 * sin(3 * w0 * t + 0.2);
vb = vpk * sin(w0 * t - 2 * pi / 3) + 5 * sin(5 * w0 * t - 0.3);
vc = vpk * sin(w0 * t + 2 * pi / 3) + 4 * sin(7 * w0 * t + 0.5);
vab = va - vb;
vbc = vb - vc;
vca = vc - va;

truth = struct('G_ab', 0.08, 'G_bc', 0.10, 'G_ca', 0.12, ...
    'model', 'delta_resistive');
iab = truth.G_ab * vab;
ibc = truth.G_bc * vbc;
ica = truth.G_ca * vca;
ia = iab - ica;
ib = ibc - iab;
ic = ica - ibc;

data = struct('t', t, 'vab', vab, 'vbc', vbc, 'vca', vca, ...
    'ia', ia, 'ib', ib, 'ic', ic, 'fs', config.fs, ...
    'f0', config.f0, 'truth', truth, ...
    'label', 'fixture three-wire delta conductance');
end
