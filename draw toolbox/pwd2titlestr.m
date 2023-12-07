function varargout = pwd2titlestr(varargin)
% titstr = pwd2titlestr(pwdstr,n)
% [titstr_1, titstr_2, ...] = pwd2titlestr(pwdstr_1,pwdstr_2,...,n)
%
% convert a directory and/or filename string for use as a string in a
% figure.
% replace '\' with '\\'
% replace '_' with '\_'
%
% n (optional, default = 0) titstr returns only the trailing part of pwdstr that 
% contains n number of '\'

nin = length(varargin);

% check for n parameter
if ~isstr(varargin{end}),
    n = varargin{end};
    nin = nin-1;
else
    n = 0;
end

for ii = 1:nin,
    titstr = varargin{ii};
    
    % trim off path
    [titpath, titname, titext] = fileparts(titstr);
    if n == 0
        titpath = {''};
    else
        titpath = split(titpath, filesep);
        n = min(n, length(titpath)-1);
        titpath = titpath(end-n:end);
    end
    titstr = fullfile(titpath{:}, [titname titext]);
    
    % search for '\' and replace with '\\'
    titstr = replace(titstr, '\', '\\');
    
    % search for '_' and replace
    titstr = replace(titstr, '_', '\_');
    
    % output
    varargout{ii} = titstr;
end

return
