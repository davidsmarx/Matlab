function [xn, x0, a] = getFHTgrid(N)
% [xn, x0, a] = getFHTgrid(N)
%
% xn = grid points
% a  = alpha parameter
% x0 = xn(1)

if nargin == 0, error('usage: [xn, x0, a] = getFHTgrid(N)'); end

a = Geta(N);
x0 = (1+exp(a)).*exp(-a*N)/2;
xn = x0.*exp(a.*[0:N-1]');  % length = N, n = 0..N-1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = Geta(N);
% solve a = -log(1-exp(-a))/(N-1)

a0 = 4/N; % initial guess
e0 = 1;
while e0 > 1e-15;
    a = -log(1-exp(-a0))/(N-1);
    e0 = abs(a-a0);
    a0 = a;
end
