function val = CheckOption(varstring, defaultval, varargin)
% val = CheckOption(varstring, defaultval, varargin)
% val = CheckOption(varstring, defaultval, struct)
%
% utility for checking varargin for an option
% varstring = string keyword
% if varstring is found in varargin{:},
% the value of the following entry in varargin is returned
% else the defaultval is returned
%
% example:
% RminSc = CheckOption('RminSc', S.RminSc, varargin{:});
val = defaultval;

if isempty(varargin)
    return
end

% options are in a struct
if isstruct(varargin{1}),
    sOpt = varargin{1};
    if isfield(sOpt, varstring),
        val = sOpt.(varstring);
    end    
end

% options are a list of name, value pairs
iv = find(strcmp(varargin, varstring));
if ~isempty(iv),
    val = varargin{iv+1};
end

