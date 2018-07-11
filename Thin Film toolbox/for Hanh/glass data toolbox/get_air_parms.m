function gout = get_air_parms(lam,temp,CTE)
% gout = get_air_parms(lam,temp)
% 
% lam = wavelength [um]
% temp = temperature [C]
% CTE = CTE of spacer material (optional)
% 
% gout is struct:
%    n is index (absolute) at lam and temp
%    cte is CTE
%    dndt is 1e6 * dndt (absolute) at lam and temp
%    strain is n*a + dndt
%    name is name

[nair, dndtair] = indexofair(lam,temp,1);

gout.n = nair(:);
if exist('CTE','var'), gout.cte = CTE; else, gout.cte = 0; end
gout.dndt = 1e6.*dndtair(:);
gout.strain = gout.n.*gout.cte + gout.dndt;
gout.name = 'air';

return