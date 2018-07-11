function y = upsample(x,a)
% y = upsample(x,a)
% a = factor (>1) by which to upsample

if nargin < 2, error('usage: y = upsample(x,a)'); end

[N Nx] = size(x);
if N == 1, x = x(:); [N Nx] = size(x); end

Ny= a*N;
nz= Ny-N;

X = myfftshift(fft(x));
Y = [zeros(floor(nz/2),Nx); X; zeros(ceil(nz/2),Nx)];
y = ifft(myfftshift(Y));

return