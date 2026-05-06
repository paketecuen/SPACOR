function value = nv3_get_option(opts, name, defaultValue)
%NV3_GET_OPTION Return struct option with a default.

if nargin < 3
    defaultValue = [];
end
if isstruct(opts) && isfield(opts, name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end

end
