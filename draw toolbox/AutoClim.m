function clim = AutoClim(A, varargin)
% clim = AutoClim(A, options)
%
% options:
%   'one-sided' [false] or true
%   'two-sided' [true] or false
%   'symmetric' [false] or true (i.e. two-sided and symmetric)
%
% A can be any dimensions, only A(:) is used
% A can be an axes handle

switch class(A)
    case 'matlab.graphics.axis.Axes',
        hax = A;
        hh = get(hax,'Children');
        hIm = hh( strcmpi(get(hh,'type'), 'image') );
        if ~isa(hIm, 'matlab.graphics.primitive.Image'),
            error('cannot interpret input A');
        end        
        A = get(hIm,'CData');
        
    case 'double'
        % do nothing
        
    otherwise
        A = double(A);
end

% options
bOnesided = CheckOption('one-sided',false,varargin{:});
bSymmetric = CheckOption('symmetric',false,varargin{:});
pctscale = CheckOption('pctscale', 100, varargin{:});

if bOnesided && bSymmetric,
    % doesn't make sense
    warning(['cannot be both one-sided and symmetric, setting to symmetric']);
    bOnesided = false;
end

% check
if isempty(A),
    error('AutoClim: input is empty');
end
if range(A(:)) <= eps*mean(abs(A(:))),
    error('AutoClim: range is 0');
end

% for complex numbers, sort sorts the abs()
asort = sort(A(:),'ascend');
amax = asort(ceil((pctscale/100)*length(asort)));
amin = asort(floor((1-(pctscale/100))*length(asort))+1);

% check result is valid
if amax - amin <= 0,
    warning(['clim = ' num2str([amin amax]) '; clim must be monotonic increasing. Switching to [min max].']);
    amin = min(A(:));
    amax = max(A(:));
end


if bOnesided,
    clim = [0 amax];

elseif bSymmetric,
    % two-sided, symmetric around a==0
    clim = max(abs([amin amax]))*[-1 1];
    
else,
    % two-sided not symmetric
    clim = [amin amax];

end
