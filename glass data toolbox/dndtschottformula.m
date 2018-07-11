function dndt = dndtschottformula(D,n,lam,dT)
% dndt = dndtschottformula(D,n,lam,dT)
% 
% D = [D0 D1 D2 E0 E1 LamTk]
% n = index at the reference temperature relative to air
% lam = wavelength [um]
% dT is temp difference from reference temperature
% result dn/dT is for absolute index
% to get dn/dT for index relative to air = (1/nair)*(dn/dT - n * dnair/dT)


dndt = 0.5*((n.^2-1)./n).*...
   (D(1) + 2*D(2).*dT + 3*D(3).*dT.^2 + (D(4) + 2*D(5).*dT)./(lam.^2-D(6).^2));

return
