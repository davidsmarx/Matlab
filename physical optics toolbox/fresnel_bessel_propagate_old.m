function [u2, rn] = fresnel_bessel_propagate_old(u1,A,lam,z,N,varargin)
% [u2, rn] = fresnel_bessel_propagate(u1,A,lam,z,N,p1,p2,...)
%
% inputs:
%    u1 can be a function handle, inline function, or vector. If u1 is a
%    vector, it must be of length N and sampled on the appropriate
%    non-uniform grid for the Fast Hankel Transform
%    p1,p2,... are extra parameters if u1 is a function, e.g.
%    u1(r,p1,p2,...)

Nf = A.^2./(lam.*z);

% test value of N, set default based on type of u1
if nargin < 5 | isempty(N), N = GetDefaultN(u1,Nf); end

% get the grid used in the FHT
[xn, x0, a] = getFHTgrid(N);
rn = A.*xn;

% calculate input function on proper grid
u1p = GetU1(u1,rn,varargin{:});

% apply quadratic phase for propagation
u1p = u1p.*exp(j.*pi.*rn.^2./(lam.*z));

% fresnel-bessel transform
[u2, xn] = FHT(u1p,Nf,N);

% apply quadratic phase correction
u2 = -j.*A.^2.*u2.*exp(j*pi.*Nf.*(xn.^2))./(lam.*z);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N = GetDefaultN(u1,Nf)

switch class(u1),
    case {'inline','function_handle'}
        N = 2.^ceil(log2(Nf)+1);
    case 'double',
        N = length(u1);
    otherwise,
        error('unrecognized class: u1');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fn = GetU1(fun,rn,varargin)

N = length(rn);

switch class(fun),
    case {'inline','function_handle'}
        fn = feval(fun,rn(1:N),varargin{:});
    case 'double',
        [nr, nc] = size(fun);
        if nr ~= N, error('input data length mismatch'); end
        fn = fun; 
    otherwise,
        error('unrecognized class: fun');
end
return
