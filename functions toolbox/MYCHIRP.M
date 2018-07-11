function y = mychirp(t,B,T)
% y = mychirp(t,B,T)
% y = mychirp(t,B,T)
% B = bandwidth (scalar)
% T = pulsewidth (scalar)

if (length(B)~=1) | (length(T)~=1),
   error('B and T must be scalars')
end

tp = ( (t<T/2) & (t>-T/2) );
y  = zeros(size(t));

y(tp) = exp(j*pi*(B/T)*t(tp).*t(tp));
