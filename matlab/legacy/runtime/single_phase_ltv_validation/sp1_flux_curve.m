function lambda = sp1_flux_curve(i, truth)
%SP1_FLUX_CURVE Saturating flux linkage lambda(i).

i = i(:);
lambda = truth.Lmin * i + (truth.L0 - truth.Lmin) * truth.Is .* ...
    tanh(i ./ truth.Is);

end
