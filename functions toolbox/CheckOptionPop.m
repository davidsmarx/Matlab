function [val, return_list] = CheckOptionPop(varstring, defaultval, varargin)
% val = CheckOptionPop(varstring, defaultval, varargin)
% val = CheckOptionPop(varstring, defaultval, struct)
%
% utility for checking varargin for an option,
% then remove the keyword, value pair from the input list
%
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
        sOpt = rmfield(sOpt, varstring);
    end
    return_sOpt = sOpt;
    list_options = varargin(2:end);
    
else
    return_sOpt = [];
    list_options = varargin;

end

% remaining options are a list of name, value pairs
iv = find(strcmp(list_options, varstring));
if ~isempty(iv),
    val = list_options{iv+1};
    % remove from list
    list_options(iv:iv+1) = [];
end

return_list = list_options;
if ~isempty(return_sOpt)
    return_list = [{return_sOpt} return_list(:)'];
end

    

