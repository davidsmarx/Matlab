function clim = AutoClim(A, varargin)
% clim = AutoClim(A, options)
%
% options:
%   'one-sided' 
%   'two-sided' (default)
%   'symmetric' (i.e. two-sided and symmetric)

% options
bOnesided = CheckOption('one-sided',false,varargin{:});
bSymmetric = CheckOption('symmetric',false,varargin{:});

if bOnesided && bSymmetric,
    % doesn't make sense
    warning(['cannot be both one-sided and symmetric, setting to symmetric']);
    bOnesided = false;
end

% for complex numbers, sort sorts the abs()
asort = sort(A(:),'ascend');
amax = asort(ceil(0.99*length(asort)));
amin = asort(floor(0.01*length(asort)));

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
