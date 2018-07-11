function y = newchirp_H(f,B,T)
% y = newchirp_H(f,B,T)
% the matched filter in frequency domain of the "newchirp"
% B = bandwidth (scalar)
% T = pulsewidth (scalar)

if (length(B)~=1) | (length(T)~=1),
   error('B and T must be scalars')
end
if B<0,
   error('B must be positive')
end

N = length(f);
df = f(2)-f(1);
dt = 1/(N*df);
t = f*(dt/df);

fc = 2*B;  % must be greater than B/2

tp = ( (t<T/2) & (t>-T/2) );
y  = zeros(size(t));

t0 = T*fc/B;
al = (2*fc+B)*(-T/2 + t0);

y(tp) = exp(j*pi*al*log(t(tp) + t0)).*exp(-j*2*pi*fc*t(tp));

y = myfftshift(fft(myfftshift(y)));
y = conj(y);

