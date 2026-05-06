function audit = run_three_wire_identifiability_audit()
%RUN_THREE_WIRE_IDENTIFIABILITY_AUDIT Check pure-sinusoid rank loss in 3-wire RL formulas.
%
% This script audits the excitation geometry behind the old 3F examples:
% - R-only delta uses a 2D area denominator and can work with one sinusoid.
% - RL/GL three-wire formulas use a 4D osculating volume; a single frequency
%   makes that volume zero, even if the two line voltages are unbalanced.

w = 2*pi*50;
t = 0:1e-5:0.04;

cases = {
    make_case("balanced_fundamental", t, w, [120, 120], [0, -2*pi/3], [], [], [])
    make_case("unbalanced_fundamental", t, w, [120, 87], [0.21, -1.83], [], [], [])
    make_case("tiny_5th_like_old_test", t, w, [120, 120], [0, -2*pi/3], 5, [0.01, 0.01], [0, -5*2*pi/3])
    make_case("akagi_like_harmonics", t, w, [sqrt(2)*130, sqrt(2)*120], [0, -2*pi/3], [7, 5], [sqrt(2)*130/100.2, sqrt(2)*120/12], [0, -2*pi/3])
    make_case("rich_old_delta_like", t, w, [120, 1.1*120], [0, -0.5*2*pi/3], [4, 16], [0.9*100, 120], [0, -4*2*pi/3])
};

rows = cell(numel(cases), 1);
for k = 1:numel(cases)
    c = cases{k};
    [s12, s1234n] = local_denominators(c, w);
    Z = [c.vab(:), c.vbc(:), c.dvab(:)/w, c.dvbc(:)/w];
    sv = svd(Z, "econ");
    tol = max(size(Z))*eps(max(sv));
    zRank = sum(sv > tol);
    rows{k} = table( ...
        c.name, ...
        zRank, ...
        rms(s12), ...
        rms(s1234n), ...
        max(abs(s1234n)), ...
        sv(end)/sv(1), ...
        'VariableNames', {'case_name','rank_v_dv','rms_s12','rms_norm_s1234','max_norm_s1234','sv_min_over_max'});
end

audit = vertcat(rows{:});
disp(audit);

outDir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
writetable(audit, fullfile(outDir, 'three_wire_identifiability_audit.csv'));
end

function c = make_case(name, t, w, baseAmp, basePhase, harmonics, harmonicAmp, harmonicPhase)
c.name = string(name);
c.vab = baseAmp(1)*cos(w*t + basePhase(1));
c.vbc = baseAmp(2)*cos(w*t + basePhase(2));
c.dvab = -w*baseAmp(1)*sin(w*t + basePhase(1));
c.dvbc = -w*baseAmp(2)*sin(w*t + basePhase(2));
c.ddvab = -w^2*baseAmp(1)*cos(w*t + basePhase(1));
c.ddvbc = -w^2*baseAmp(2)*cos(w*t + basePhase(2));
c.dddvab = w^3*baseAmp(1)*sin(w*t + basePhase(1));
c.dddvbc = w^3*baseAmp(2)*sin(w*t + basePhase(2));
c.ddddvab = w^4*baseAmp(1)*cos(w*t + basePhase(1));
c.ddddvbc = w^4*baseAmp(2)*cos(w*t + basePhase(2));

if isempty(harmonics)
    return;
end

if isscalar(harmonics)
    harmonics = [harmonics, harmonics];
end
if isscalar(harmonicAmp)
    harmonicAmp = [harmonicAmp, harmonicAmp];
end
if isempty(harmonicPhase)
    harmonicPhase = zeros(1, numel(harmonics));
end
if isscalar(harmonicPhase)
    harmonicPhase = [harmonicPhase, harmonicPhase];
end

for idx = 1:numel(harmonics)
    h = harmonics(idx);
    if idx == 1
        c.vab = c.vab + harmonicAmp(idx)*cos(h*w*t + harmonicPhase(idx));
        c.dvab = c.dvab - h*w*harmonicAmp(idx)*sin(h*w*t + harmonicPhase(idx));
        c.ddvab = c.ddvab - (h*w)^2*harmonicAmp(idx)*cos(h*w*t + harmonicPhase(idx));
        c.dddvab = c.dddvab + (h*w)^3*harmonicAmp(idx)*sin(h*w*t + harmonicPhase(idx));
        c.ddddvab = c.ddddvab + (h*w)^4*harmonicAmp(idx)*cos(h*w*t + harmonicPhase(idx));
    else
        c.vbc = c.vbc + harmonicAmp(idx)*cos(h*w*t + harmonicPhase(idx));
        c.dvbc = c.dvbc - h*w*harmonicAmp(idx)*sin(h*w*t + harmonicPhase(idx));
        c.ddvbc = c.ddvbc - (h*w)^2*harmonicAmp(idx)*cos(h*w*t + harmonicPhase(idx));
        c.dddvbc = c.dddvbc + (h*w)^3*harmonicAmp(idx)*sin(h*w*t + harmonicPhase(idx));
        c.ddddvbc = c.ddddvbc + (h*w)^4*harmonicAmp(idx)*cos(h*w*t + harmonicPhase(idx));
    end
end
end

function [s12, s1234n] = local_denominators(c, w)
s12 = c.vab.*c.dvbc - c.vbc.*c.dvab;
n = numel(c.vab);
s1234n = zeros(n, 1);
for idx = 1:n
    M = [
        c.vab(idx),        c.vbc(idx),        c.dvab(idx)/w,       c.dvbc(idx)/w
        c.dvab(idx)/w,     c.dvbc(idx)/w,     c.ddvab(idx)/w^2,    c.ddvbc(idx)/w^2
        c.ddvab(idx)/w^2,  c.ddvbc(idx)/w^2,  c.dddvab(idx)/w^3,   c.dddvbc(idx)/w^3
        c.dddvab(idx)/w^3, c.dddvbc(idx)/w^3, c.ddddvab(idx)/w^4,  c.ddddvbc(idx)/w^4
    ];
    rowScale = prod(vecnorm(M, 2, 2));
    s1234n(idx) = det(M) / max(rowScale, eps);
end
end
