function [rb, ra, tb, ta] = TFfilter2z(n,theta,tetm)
% [rb, ra, tb, ta] = TFfilter2z(n,theta,tetm)
% n = index of refraction for each layer. 
%     n(1) = index of substrate
%     n(N) = index of incident medium
%     then length(n) must be >= 2
% theta = angle of incidence, scalar, [rad]
% tetm: 0 => TE (default), 1 => TM
%ra = denominator of transfer function
%rb = numerator of transfer function
%ta = denominator of transfer function
%tb = numerator of transmission response

N = length(n);
%n = n(:);

kx = n(1)*sin(theta);
kz = -sqrt( n.^2 - kx.^2 ); % sign agrees with measurement convention
if tetm == 1,
   kzz = kz./(n.^2);
else,
   kzz = kz;
end

i0 = 1:N-1;
i1 = 2:N;
rf  = (kzz(i1) - kzz(i0))./(kzz(i1) + kzz(i0));
tf  = 2*kzz(i0)./(kzz(i1) + kzz(i0));

[rb,ra,tb,ta] = rcoeffs2tf(rf);


return
