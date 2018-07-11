R = 300*MM; % radius of curvature
A = 10*MM;  % radius of model
a = 5*MM;   % radius of beam

N = 2^19;   % number of grid points

lam = 1.3*UM; % wavelength
k   = 2*pi./lam;

rplotmax = 10*lam/(a/(R/2));
rn       = A*getFHTgrid(N);
nplot    = rn<=rplotmax;
rplot    = [-flipud(rn(nplot)); rn(nplot)]/MM;

% startbeam = inline('rect(r/(2*a)).*exp(-j*k*sqrt( (R./2 - 0.25*r.^2./R).^2 + r.^2 ))','r','R','a','k');
% idealbeam = inline('rect(r/(2*a)).*exp(-j*k*sqrt( (R./2).^2 + r.^2 ))','r','R','a','k');
startbeam = inline('rect(r/(2*a)).*exp(-j*k.*(r.^2./R + 8.4*(2*pi/k).*(r./a).^4))','r','R','a','k');
idealbeam = inline('rect(r/(2*a)).*exp(-j*k.* r.^2./R                           )','r','R','a','k');

zlist = linspace(R/2-1*MM,R/2,11)';

for ii = 1:length(zlist),
    z = zlist(ii);
    
    [u1, rn] = fresnel_bessel_propagate(startbeam,A,lam,z,N,R,a,k);
%    [u1, rn] = fresnel_bessel_propagate(idealbeam,A,lam,z,N,R,a,k);
    
%     atmp = interp1(rn,abs(u1),rr);
    atmp = abs(u1(nplot));

    magu{ii} = [flipud(atmp); atmp];

end

figure, plot(rplot,[magu{:}]), grid
