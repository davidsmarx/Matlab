function U = blackbodyspectrum(lam,tt)
% U = blackbodyspectrum(lam,T)
%
% Spectral energy density, U(?,T) = 8?hc?^-5 /  ( e^(hc/?kT)-1 ) 
%
% lam = wavelength, (m)
% T is black body temperature (Kelvin)



global C HPLANCK KBOLTZMANN;

U1 = 8.*pi.*HPLANCK.*C.*(lam.^-5);
U2 = exp( HPLANCK.*C./(lam.*KBOLTZMANN.*tt) ) - 1;
U = U1 ./ U2;
