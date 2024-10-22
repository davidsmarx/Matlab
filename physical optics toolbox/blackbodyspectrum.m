function U = blackbodyspectrum(lam,tt)
% U = blackbodyspectrum(lam,T)
%
% Spectral energy density, U(lam,T) = 8 pi h c lam^-5 /  ( e^(h c/lam k T) - 1 ) 
%
% lam = wavelength, (m)
% T is black body temperature (Kelvin)



global C HPLANCK KBOLTZMANN;

U1 = 8.*pi.*HPLANCK.*C.*(lam.^-5);
U2 = exp( HPLANCK.*C./(lam.*KBOLTZMANN.*tt) ) - 1;
U = U1 ./ U2;
