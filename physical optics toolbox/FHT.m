function [g, xn] = FHT(fun,Nf,N,varargin)
% [g, xn] = FHT(fun,Nf,N,p1,p2,...)
% 
% Fast Hankel Transform of fun(x,p1,p2,...)
%
%  FHT(y) = 2 pi int( fun(x) J0(2 pi Nf y x) x dx, [0 1] )
%
% input fun: function handle, inline function, or vector of samples sampled
% on the proper xn grid. For example a vector returned by a previous call
% to FHT.
% Nf = fresnel number, = A^2/(lam*z), where A is the upper integration
% limit, lam = wavelength, z = propagation distance.
% N = size of the transform. Choose N > Nf, and N = 2^m, for integer m.
%
% outputs:
%   g = Hankel transform of fun
%   xn = sample points where g is evaluated. The sample grid is not uniform
%
% algorithm:
%   Magni, Cerullo, and De Silvestri, "High-accuracy Fast Hankel Transform
%   for Optical Beam Propagation," JOSA-A, Vol. 9, No. 11, p. 2031.

% non-uniform sampling
[xn, x0, a] = getFHTgrid(N); % length = N, n = 0..N-1

% evaluate the function (or validate input vector)
fn = Getfn(fun,xn,varargin);

% pn: length 2N, p(n), n=0..2N-1
pn = [-diff(fn).*exp(a.*([0:N-1]'+1-N)); zeros(N,1)];
k0 = (2.*exp(a) + exp(2*a))./( (1+exp(a)).^2.*(1-exp(-2*a)) );
pn(1) = k0.*pn(1);

% jn: length 2N, n=0..2N-1
jn = besselj(1,2.*pi.*Nf.*x0.*exp(a.*([0:2*N-1]'+1-N)));

%g = fft(fft(pn).*ifft(jn));
pn = fft(pn);
jn = ifft(jn);
pn = jn.*pn;
g = fft(pn);

g = g(1:N)./(Nf.*xn);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fn = Getfn(fun,xn,funparms)

N = length(xn);

switch class(fun),
    case {'inline','function_handle'}
        fn = [feval(fun,xn(1:N),funparms{:}); 0]; % length = N+1, f(xn), n=0..N
    case 'double',
        [nr, nc] = size(fun);
        if nr ~= N, error('input data length mismatch'); end
        fn = [fun; 0]; % length = N+1
    otherwise,
        error(['unrecognized class fun: ' class(fun)]);
end
return
