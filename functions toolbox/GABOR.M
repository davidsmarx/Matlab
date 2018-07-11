function [G, xsh] = gabor(p,h,Nx)
% G = gabor(p,h,Nx)
% p = real data vector
% h = window, length(h) must be < length(p)
% Nx= number of time steps to use
% xsh = list of time steps that were used

if ~isreal(p), error('input data vector must be real'); end

Nh = length(h);
Np = length(p);
if Np < Nh, error('length of window h must be less than length of data'); end

Nfft = 2.^(ceil(log2(2*Nh-1)));
fprintf('Nfft = %d\n',Nfft);

xsh= round(0:(Np/Nx):Np-1);
if length(xsh)~=Nx fprintf('length(xsh)=%d Nx=%d\n',length(xsh),Nx); error('xshift'); end

Nx = length(xsh);
G  = zeros([Nx Nfft]);

py = zeros([Np+Nfft 1]); py(1:Np) = p; % pad end of input data with Nfft zeros
for i = 1:Nx,
   x  = xsh(i);
   pp = shift(py,x);
   z  = pp(1:Nh);
   Z  = ifft(h .* (z-mean(z)),Nfft);
   G(i,:) = Z(1:Nfft)'; % only need to take the positive frequencies since p is real
end

%figure, waterfall(([-Nfft/2:Nfft/2-1])/Nfft,xsh,abs(myfftshift(G'))')

return

function y = shift(x,m)

n = length(x);
y = zeros(size(x));
y(1:n-m)   = x(m+1:n);
y(n-m+1:n) = x(1:m);
