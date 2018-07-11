function ee = laguerre_gauss_beam(r,t,z,p,l,w0,lam)
% ee = laguerre_gauss_beam(r,p,z,p,l,w0,lam)
%
% r,t,z are cylindrical coordinates
%
% p = radial mode number
% l = azimuthal mode number

k  = 2*pi./lam;
zr = pi.*w0.^2./lam;
w2 = w0.^2.*(1 + (z./zr).^2);
if abs(z) < 1e-6*lam,
    rr = inf;
else,
    rr  = z.*(1 + (zr./z).^2);
end


Lpl = mfun('orthopoly[L]',p,l,2.*r.^2./w2);
% or
Lpoly_pl = inline(maple(['orthopoly[L](' num2str(p) ',' num2str(l) ',x)']),'x')
% Lpoly_pl(2.3)

amp = (sqrt(2/w2).*r).^l .* Lpl .* exp(-r.^2./w2);
phase = k*z - (2*p+l+1).*atan(z./zr) + 0.5*k.*r.^2./rr + l.*t;

ee = amp.*exp(j*phase);