function U = FresnelDiffRect(x,y,lx,ly,lam,z)
% U = FresnelDiffRect(x,y,lx,ly,lam,z)
%
% x,y = vector of x,y coordinates to evaluate the field
% lx, ly = width, height of rectangle, centered at the origin
% lam = wavelength
% z = propagation distance
% all values are the same unit (meters)
%
% U = matrix of complex field values at (X,Y) (the geometric phase,
% exp(jkz) is not included)
%
% formulation from Goodman, eqns. 4-27 to 4-30

K = 2*pi/lam;
aa = sqrt(K./pi./z);

zet1 = -aa.*(lx/2 + x);
zet2 = aa.*(lx/2 - x);

eta1 = -aa.*(ly/2 + y);
eta2 = aa.*(ly/2 - y);

% define FresnelC and FresnelS as inline functions for ease
fc = inline('mfun(''FresnelC'',x)','x');
fs = inline('mfun(''FresnelS'',x)','x');

Ux = fc(zet1)-fc(zet2) + j*( fs(zet1)-fs(zet2) );
Uy = fc(eta1)-fc(eta2) + j*( fs(eta1)-fs(eta2) );

if length(x) > 1 & length(y) > 1,
    % calculate the field on a 2-d grid of x,y points
    U = -(0.5*j) .* (Uy(:) * Ux(:).');
else
    U = -(0.5*j) .* Uy .* Ux;
end
