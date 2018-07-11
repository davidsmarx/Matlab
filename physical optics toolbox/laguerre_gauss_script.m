unitsdefinitions;

w0  = 0.1*MM;
z   = 1*UM;
lam = 1.0*UM;

r = linspace(0,4*w0);
phi = linspace(-pi,pi);
[R, Phi] = meshgrid(r,phi);

zlist = [-3:3]*15*MM;

for ii = 1:length(zlist),
    z = zlist(ii);
p = 1;
l = 2;
e1 = laguerre_gauss_beam(R,Phi,z,p,l,w0,lam);

p = 1;
l = 1;
e2 = laguerre_gauss_beam(R,Phi,z,p,l,w0,lam);

ee = e1 + e2;

[X, Y] = deal(R.*cos(Phi),R.*sin(Phi));
figure, surf(X/MM,Y/MM,abs(ee)),...
    set(gca,'view',[90,90]), set(gca,'dataaspectratio',[1 1 2.5])
% figure, surf(X/MM,Y/MM,angle(ee)/pi),...
%     set(gca,'viennnw',[90,90]), set(gca,'dataaspectratio',[1 1 2.5])

end