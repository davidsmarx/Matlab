function [opd, pow, phasor] = Heterodyne(varargin)
% [opd, power, phasor] = Heterodyne(Usig,Uref,rn,D,lam)
% [opd, power, phasor] = ...
%     Heterodyne(sigbeam, refbeam, detshape, detD, x, y, dx, dy, lam)
%
% D = diameter for detector integration

switch length(varargin),
    case 5,
        [opd, pow, phasor] = Heterodyne_FHT(varargin{:});
    case 9,
        [opd, pow, phasor] = Heterodyne_BB(varargin{:});
    otherwise,
        error('error calling Heterodyne');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opd, power, phasor] = Heterodyne_BB(varargin)

[sigbeam, refbeam, detshape, detD, x, y, dx, dy, lam] = deal(varargin{:});

K = 2*pi/lam;

hetamp = sigbeam .* conj(refbeam);

[X, Y] = meshgrid(x,y);
switch detshape,
    case 'circle',
        R = sqrt(X.^2 + Y.^2);
        hetamp(R>detD/2) = 0;
    case 'square',
        hetamp(Y>detD/2&Y<-detD/2) = 0;
        hetamp(X>detD/2&X<-detD/2) = 0;
    otherwise,
        error(['detector shape ' detshape ' not supported']);
end
phasor = sum(hetamp(:))*dx*dy;

power = abs(phasor);
opd = angle(phasor)./K;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opd, pow, Eh] = Heterodyne_FHT(varargin)

[Usig, Uref, rn, D, lam] = deal(varargin{:});

K = 2*pi./lam;

nd = rn<(D/2);

Eh = 2*pi*sum( Usig(nd).*conj(Uref(nd)) .*rn(nd).*gradient(rn(nd)));
opd = angle(Eh)./K;
pow = abs(Eh);

return