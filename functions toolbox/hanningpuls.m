function s = hanningpuls(t,T)
% s = hanningpuls(t,T)
% return a pulse formed from a Hanning window with width T
% t = time is a real 1-d vector

tp = (t>=-T/2)&(t<T/2);
s = zeros(size(t));

s(tp) = 0.5*(1.0 + cos(2*pi*t(tp)/T));