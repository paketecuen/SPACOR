function Linc = sp1_incremental_inductance(i, truth)
%SP1_INCREMENTAL_INDUCTANCE Incremental inductance d lambda / d i.

i = i(:);
u = i ./ truth.Is;
Linc = truth.Lmin + (truth.L0 - truth.Lmin) ./ cosh(u) .^ 2;

end
