function [u2, rn] = fresnel_bessel_propagate(u1,A,lam,z,N,varargin)
% [u2, rn] = fresnel_bessel_propagate(u1,A,lam,z,N,p1,p2,...)
%
% inputs:
%    u1 can be a function handle, inline function, or vector. If u1 is a
%    vector, it must be sampled on the appropriate non-uniform grid for 
%    the Fast Hankel Transform.
%    A = maximum extent of sampling grid. If u1 is a vector, A = radius for
%    u1(end).
%    N is the length of the sampling grid. If u1 is a vector,
%    then N is set to length(u1), regardless of input. If u1 is a function,
%    and N = [], then N is chosen based on the fresnel number.
%    p1,p2,... are extra parameters if u1 is a function, e.g.
%    u1(r,p1,p2,...)

[u1p,xn,rn,N,Nf] = ParseInputs(u1,A,lam,z,N,varargin{:});
% figure, plot(A*xn,abs(u1p),'-o'), grid, title('source beam amplitude')
%disp(['Nf = ' num2str(Nf)]);

% apply quadratic phase for propagation
u1p = u1p.*exp(j.*pi.*rn.^2./(lam.*z));

% fresnel-bessel transform
[u2, xn] = FHT(u1p,Nf,N);

% apply quadratic phase correction
u2 = -j.*A.^2.*u2.*exp(j*pi.*Nf.*(xn.^2))./(lam.*z);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u1p,xn,rn,N,Nf] = ParseInputs(u1,A,lam,z,Nin,varargin)

% test valid A input
if isempty(A) | A<=0, error('A must be > 0 for functional u1'); end
% calculate fresnel number
Nf = A.^2./(lam.*z);

switch class(u1),
    case {'inline','function_handle'}
        % test valid N input
        if isempty(Nin),
            N = 2.^ceil(log2(Nf)+1);
        else,
            N = Nin;
        end
        
        % get the grid used in the FHT
        [xn, x0, a] = getFHTgrid(N);
        rn = A.*xn;

        % evaluate function
        u1p = feval(u1,rn,varargin{:});
        
    case 'double',
        [nr, nc] = size(u1);
        N = nr;

        % get the grid used in the FHT
        [xn, x0, a] = getFHTgrid(N);
        rn = A.*xn;

        u1p = u1; 
    otherwise,
        error('unrecognized class: u1');
end
