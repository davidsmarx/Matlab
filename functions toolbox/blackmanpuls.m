function s = blackmanpuls(t,T)
% s = blackmanpuls(t,T)
% return a pulse formed from a Blackman window with width T
% t = time is a real 1-d vector

tp = (t>=-T/2)&(t<T/2);
s = zeros(size(t));

s(tp) = 0.42 + 0.50*cos(2*pi*t(tp)/T) + 0.08*cos(4*pi*t(tp)/T);