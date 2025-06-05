function gout = get_hpfs_parms(lam,temp)
% gout = get_hpfs_parms(lam,temp)
% lam [m]
% temp [deg C]
%
% gout is a struct:
%   n, cte, dndt, strain, name

% Reference Conditions
Tref = 22; % [C]
% index measured in 760mm Hg Dry Nitrogen Gas

% HPFS constants (from Corning® HPFS® 7979, 7980, 8655 Fused Silica
%    Optical Materials Product Information)

% Sellmeier Dispersion Equation HPFS 7980 @ 22 C
P = [0.683740494 0.00460352869 0.420323613 0.0133968856 0.58502748 64.4932732];

% dn/dT dispersion equation (not the Schott formula)
D = [9.39059 0.23529 -1.31856e-3 3.02887e-4];

% CTE for 5 to 35 C
cte = 0.52; % 10^-6

% index relative to Dry Nitrogen Gas
n = dispersionformula2index(P,'Sellmeier1',lam);

% dn/dT equation is from Corning reference, lam in um
% dn/dT is [ppm/C]
lam2 = (lam./1e-6).^-2;
dndt = D(1) + lam2.*(D(2) + lam2.*(D(3) + lam2.*D(4))); % [1e-6]

% index relative to air at temp
nabs = n + 1e-6.*(temp-Tref).*dndt;
nair = indexofair(lam,temp,1);

gout.n = nabs; %nabs(:)./nair(:);
gout.cte = cte;
gout.dndt = dndt;
gout.strain = gout.n.*gout.cte + gout.dndt;
gout.name = 'HPFS';

return