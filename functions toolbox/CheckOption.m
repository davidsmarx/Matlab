function val = CheckOption(varstring, defaultval, varargin)
% CheckOption(varstring, defaultval, varargin)
%
% utility for checking varargin for an option
% if found, the following entry in varargin is returned
% else the defaultval is returned

iv = find(strcmp(varargin, varstring));

if ~isempty(iv),
    val = varargin{iv+1};
else
    val = defaultval;
end
