function [TAchC,LAchC,TchC] = ParaxChromAberrations(obj_h,obj_d,pup_h,t,c,n1,n2)
% [TAchC,LAchC,TchC] = ParaxChromAberrations(obj_h,obj_d,pup_h,t,c,n1,n2)
%
% obj_h = object height
% obj_d = distance from object to first surface
% pup_h = pupil height
% t  = list of thicknesses, length(d) = # surfaces - 1
% c  = list of surface curvatures, length(c) = # surfaces
% n1 = list of medium indexes for the short wavelength, length(n) = # surfaces + 1
% n2 = list of medium indexes for the long wavelength, length(n) = # surfaces + 1
% 
% TAchC = paraxial transverse axial chromatic aberration
% LAchC = paraxial longitudinal axial chromatic aberration
% TchC = paraxial lateral chromatic aberration

% ray trace at mean wavelength
n = mean(n1(:)',n2(:)');

% dispersion = n1 - n2
dn = n1 - n2;

% axial ray
%u0a = p1./l1; y0a = p1;
aa = [pup_h/obj_d; pup_h];
[ua,ya,iia] = paraxraytrace_surf(t,c,n,aa);
li = -ya(end)/ua(end);

% principle ray
%u0p = h0./l1; y0p = 0;
ap = [obj_h/obj_d; 0];
[up,yp,iip] = paraxraytrace_surf(t,c,n,ap);
hi = yp(end) + li*up(end);

invar = ap(2)*n(1)*aa(1) - aa(2)*n(1)*ap(1);
hcheck = invar./(n(end)*ua(end));

nn = n(1:end-1); np = n(2:end);
dnn = dn(1:end-1); dnp = dn(2:end);

Ba = nn.*(np-nn).*ya.*(ua + iia)./(2.*np.*invar);
Bp = nn.*(np-nn).*yp.*(up + iip)./(2.*np.*invar);

TAchC = -ya.*iia.*(dnn - nn.*dnp./np)./(np(end).*ua(end));
TchC = -ya.*iip.*(dnn - nn.*dnp./np)./(np(end).*ua(end));
LAchC = -TAchC./ua(end);

return
