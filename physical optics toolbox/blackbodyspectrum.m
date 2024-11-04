function U = blackbodyspectrum(lam,tt)
% U = blackbodyspectrum(lam,T)
%
% Spectral energy density [J/m^4], U(lam,T) = 8 pi h c lam^-5 /  ( e^(h c/(lam kT)) - 1 ) 
% Spectral radiance [W / m^2 / sr / Hz] = 2 h f lam^-2 / ( e^(h c/(lam kT)) - 1)
% Spectral radiance [W / m^2 / sr / Hz] = 2 h c lam^-3 / ( e^(h c/(lam kT)) - 1)
%
% lam = wavelength, (m)
% T is black body temperature (Kelvin)



global C HPLANCK KBOLTZMANN;

% spectral energy density:
% U1 = 8.*pi.*HPLANCK.*C.*(lam.^-5);

% spectral radiance:
U1 = 2.*HPLANCK.*C.*(lam.^-3);

U2 = exp( HPLANCK.*C./(lam.*KBOLTZMANN.*tt) ) - 1;
U = U1 ./ U2;
