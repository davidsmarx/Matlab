function gout = get_hpfs_parms(lam,temp)
% gout = get_hpfs_parms(lam,temp)
%
% gout is a struct:
%   n, cte, dndt, strain, name

% Reference Conditions
Tref = 22; % [C]
% index measured in 760mm Hg Dry Nitrogen Gas

% HPFS constants:
% Sellmeier Dispersion Equation
P = [0.683740494 0.00460352869 0.420323613 0.0133968856 0.58502748 64.4932732];

% dn/dT dispersion equation (not the Schott formula)
D = [9.39059 0.23529 -1.31856e-3 3.02887e-4];

% CTE for 5 to 35 C
cte = 0.52; % 10^-6

% index relative to Dry Nitrogen Gas
n = indexformula(P,'Sellmeier1',lam);

% 
l2 = lam.^-2;
dndt = D(1) + l2.*(D(2) + l2.*(D(3) + l2.*D(4))); % [1e-6]

% index relative to air at temp
nabs = n + 1e-6.*(temp-Tref).*dndt;
nair = indexofair(lam,temp,1);

gout.n = nabs; %nabs(:)./nair(:);
gout.cte = cte;
gout.dndt = dndt;
gout.strain = gout.n.*gout.cte + gout.dndt;
gout.name = 'HPFS';

return