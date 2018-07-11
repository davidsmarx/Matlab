function varargout = pwd2titlestr(varargin)
% titstr = pwd2titlestr(pwdstr,n)
% [titstr_1, titstr_2, ...] = pwd2titlestr(pwdstr_1,pwdstr_2,...,n)
%
% convert a directory and/or filename string for use as a string in a
% figure.
% replace '\' with '\\'
% replace '_' with ' '
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
    %ns = findstr(titstr,'\');
    ns = regexp(titstr,'[\\/]');
    if ~isempty(ns),
        if n < 0,
            error('n must be >= 0');
        elseif n >= length(ns),
            % use the whole thing, do nothing
        else % n >= 0,
            titstr = titstr(ns(end-n)+1:end);
        end
    end

    % search for '\' and replace with '\\'
    ns = findstr(titstr,'\');
    if ~isempty(ns),
        tmpstr = [];
        ic = 1;
        for in = 1:length(ns),
            tmpstr = [tmpstr titstr(ic:ns(in)) '\'];
            ic = ns(in) + 1;
        end
        tmpstr = [tmpstr titstr(ic:end)];
        
        titstr = tmpstr;
    end

    
    % search for '_' and replace
    nu = findstr(titstr,'_');
    titstr(nu) = ' ';
    
    % output
    varargout{ii} = titstr;
end

return
