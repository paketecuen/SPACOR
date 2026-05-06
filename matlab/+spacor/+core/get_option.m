function value = get_option(opts, name, defaultValue)
%GET_OPTION Read an option field with a default value.

if nargin < 3
    defaultValue = [];
end
if isstruct(opts) && isfield(opts, name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end
end

