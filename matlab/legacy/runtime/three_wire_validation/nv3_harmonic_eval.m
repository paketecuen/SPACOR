function out = nv3_harmonic_eval(components, t, f0, maxOrder)
%NV3_HARMONIC_EVAL Evaluate harmonic signal, derivatives, and primitive.
%
% components rows are [amplitude harmonic phase].
% out.values{1} is the signal, out.values{k+1} is derivative order k.
% out.primitive is the zero-mean harmonic primitive.

if nargin < 4
    maxOrder = 1;
end

t = t(:);
omega0 = 2 * pi * f0;
values = cell(maxOrder + 1, 1);
for order = 0:maxOrder
    values{order + 1} = zeros(size(t));
end
primitive = zeros(size(t));

for row = 1:size(components, 1)
    amplitude = components(row, 1);
    harmonic = components(row, 2);
    phase = components(row, 3);
    omega = omega0 * harmonic;
    angle = omega * t + phase;
    for order = 0:maxOrder
        values{order + 1} = values{order + 1} + ...
            amplitude * omega ^ order * cos(angle + order * pi / 2);
    end
    primitive = primitive + amplitude / omega * sin(angle);
end

out = struct();
out.values = values;
out.primitive = primitive;

end
