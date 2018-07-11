function [rb, ra, ta] = TFfilter2z_analytic(n,tetm)
% [R, T, rr, tt] = thin_film_filter_2(n,d,theta,lam,[tetm])
% n = index of refraction for each layer. 
%     n(1) = index of incident medium
%     n(N) = index of transmission medium
%     then length(n) must be >= 2
% d = thickness of each layer, not counting incident medium or transmission
%     medium. length(d) = length(n)-2
% theta = angle of incidence [rad]
% lam = wavelength. units of lam must be same as d
% tetm: 0 => TE (default), 1 => TM
% R = normalized reflected intensity coefficient
% T =        "   transmitted     "
% rr = complex field reflection coefficient
% tt =      "        transmission "

N = length(n);

try tetm, catch tetm = 0; end

syms z;

kz = -n; % sign agrees with measurement convention

if tetm == 1,
   kzz = kz./(n.^2);
else,
   kzz = kz;
end

i1 = 1:N-1;
i2 = 2:N;
T  = 0.5*(kzz(i1) + kzz(i2))./kzz(i1);
R  = (kzz(i1) - kzz(i2))./(kzz(i1) + kzz(i2));

A  = eye(2);
for i = 1:N-1,
	Ai = T(i)*( [1 R(i);R(i)*z z] * A );
	A  = Ai;
	
end

ra = collect(A(2,2),z);   % denominator of transfer function
rb = collect(A(1,2),z);   % numerator of transfer function
ta = collect(A(2,2),z);   % denominator of transfer function

return
