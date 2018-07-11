function gout = get_schottglass_parms(ctmp,lam,temp)
% gout = get_schottglass_parms(ctmp,lam,temp)
% 
% ctmp = struct returned from read_schott_cat
% lam = wavelength [um]
% temp = temperature [C] (optional, default = 23)
% 
% gout is struct:
%    n is ABSOLUTE index at lam and temp
%    cte is CTE
%    dndt is 1e6 * dndt (absolute) at lam and temp
%    strain is n*a + dndt
%    name is name

if ~exist('temp','var'), temp = 23; end

[nair, dndtair] = indexofair(lam,temp,1);

gout.n = dispersionformula2index(...
   ctmp.disp_poly,'Sellmeier1',lam,ctmp.dndt_poly,temp).*nair(:);

if ~isempty(ctmp.alpha30_70),
   gout.cte = ctmp.alpha30_70.*ones(size(gout.n));
else,
   gout.cte = 0;
end

if any(ctmp.dndt_poly),
   gout.dndt = 1e6*dndtschottformula(ctmp.dndt_poly,gout.n./nair(:),lam,temp-20);
else,
   gout.dndt = 0;
end

gout.strain = gout.n.*gout.cte + gout.dndt;

gout.name = ctmp.name;
return