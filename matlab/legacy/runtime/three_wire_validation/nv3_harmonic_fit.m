function fit = nv3_harmonic_fit(x, t, f0, harmonics, opts)
%NV3_HARMONIC_FIT Fit x(t) with DC + harmonic basis.
%
% The primitive intentionally ignores the DC term. In electrical data that
% term is usually measurement offset; integrating it would create a ramp.

if nargin < 5
    opts = struct();
end

x = x(:);
t = t(:);
harmonics = harmonics(:).';
omega0 = 2 * pi * f0;

B = ones(numel(t), 1);
for h = harmonics
    B = [B, cos(h * omega0 * t), sin(h * omega0 * t)]; %#ok<AGROW>
end

s = svd(B, 'econ');
if isempty(s) || min(s) <= eps(max(s))
    basisCondition = Inf;
else
    basisCondition = max(s) / min(s);
end

coef = B \ x;
maxOrder = nv3_get_option(opts, 'maxOrder', 1);
values = cell(maxOrder + 1, 1);
for order = 0:maxOrder
    values{order + 1} = zeros(size(t));
end
values{1} = values{1} + coef(1);
primitive = zeros(size(t));

col = 2;
for h = harmonics
    a = coef(col);
    b = coef(col + 1);
    omega = h * omega0;
    angle = omega * t;
    for order = 0:maxOrder
        if order == 0
            term = a * cos(angle) + b * sin(angle);
        elseif mod(order, 4) == 1
            term = -a * omega ^ order * sin(angle) + b * omega ^ order * cos(angle);
        elseif mod(order, 4) == 2
            term = -a * omega ^ order * cos(angle) - b * omega ^ order * sin(angle);
        elseif mod(order, 4) == 3
            term = a * omega ^ order * sin(angle) - b * omega ^ order * cos(angle);
        else
            term = a * omega ^ order * cos(angle) + b * omega ^ order * sin(angle);
        end
        values{order + 1} = values{order + 1} + term;
    end
    primitive = primitive + a / omega * sin(angle) - b / omega * cos(angle);
    col = col + 2;
end

fit = struct();
fit.values = values;
fit.primitive = primitive;
fit.coefficients = coef;
fit.basisCondition = basisCondition;
fit.harmonics = harmonics;
fit.reconstructionResidual = norm(x - values{1}) / max(norm(x), eps);

end
