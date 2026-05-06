function value = sp1_robust_cv(x)
%SP1_ROBUST_CV Robust coefficient of variation based on IQR/median.

x = x(isfinite(x));
if isempty(x)
    value = Inf;
else
    value = iqr_local(x) / max(abs(median(x)), eps);
end

end

function value = iqr_local(x)
x = sort(x(:));
value = percentile_local(x, 75) - percentile_local(x, 25);
end

function value = percentile_local(x, p)
if isempty(x)
    value = NaN;
    return;
end
pos = 1 + (numel(x) - 1) * p / 100;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    value = x(lo);
else
    value = x(lo) + (x(hi) - x(lo)) * (pos - lo);
end
end
