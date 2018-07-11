unitsdefinitions;

w0  = 0.1*MM;
z   = 1*UM;
lam = 1.0*UM;

k   = 2*pi./lam;
zr  = pi.*w0.^2./lam;
w2  = w0.^2.*(1 + (z./zr).^2);
rr  = z.*(1 + (zr./z).^2);

E00 = inline('exp(-r.^2./w2).*exp(-j*(k.*z - atan(z./zr) + 0.5*k.*r.^2./rr))',...
    'r','z','w2','rr','zr','k');

Epl = inline('

r = linspace(-1*MM,1*MM,1024)';
plotampphase(r/MM,E00(r,z,w2,rr,zr,k));

p = 3;
l = 1;

mfun('orthopoly[L]',p,l,2.3)


Lpoly_pl = inline(maple(['orthopoly[L](' num2str(p) ',' num2str(l) ',x)']),'x')
Lpoly_pl(2.3)