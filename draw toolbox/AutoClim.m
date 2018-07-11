function clim = AutoClim(A, varargin)
% clim = AutoClim(A, options)
%
% options:
%   'one-sided' (default)
%   'two-sided'
%   'symmetric' (i.e. two-sided and symmetric)

% options

% for complex numbers, sort sorts the abs()
asort = sort(A(:),'ascend');
amax = asort(ceil(0.99*length(asort)));
amin = asort(floor(0.01*length(asort)));

% default is one-sided
clim = [0 amax];

if any(ismember(varargin, 'one-sided')),
    clim = [0 amax];

elseif any(ismember(varargin, 'two-sided')) ...
        || any(ismember(varargin, 'symmetric')),
    clim = [amin amax];

end

if any(ismember(varargin, 'symmetric')),
    clim = max(abs(clim))*[-1 1];    
end

