function [eo, ee, es, ep, to, te, rs, rp] = uniaxial_match_iso2uni(k1,ko,ke,se,n,m,c,ei)
% [eo, ee, er, to, te, r] = uniaxial_refraction(k1,ko,ke,se,n,kr,ei,er)
% k1 direction vector of incident ray
% ko, ke, se direction vectors ord. ray, ext. phase, ext. ray
% n direction vector of surface normal
% m column vector of indexes: [n_incident_medium n_ordinary n_extra-ordinay]';
% c crystal axis direction vector
% ei direction vector of incident polarization
% outputs:
%  eo, ee, er E-field direction vectors
%  to, te transmission coefficients
%  r reflection coefficient
% resulting E-fields are:
%   Eo = sqrt(m(2))*to*eo;
%   Ee = sqrt(m(3))*te*ee;
%   Er = sqrt(m(1))*r*er;
n1 = m(1); no = m(2); ne = m(3);

kr = k1 - 2*(n'*k1)*n; % propagation direction of reflected field
ts = cross(k1,n); if ts'*ts < 1.0e-12, ts = [0 1 0;0 0 1;1 0 0]*n; else ts = ts./sqrt(ts'*ts); end
tp = cross(n,ts); if abs( (tp'*tp)-1.0 ) > 1.0e-12, error('something wrong with tp'); end

ei = unitize(ei);
es = ts;
es = cross(kr,n); if es'*es < 1.0e-12, es = -ts; else es = es./sqrt(es'*es); end
ep = unitize(cross(kr,es));
ee = unitize(cross(se, cross(se,ke))); 
if ee'*ee < 1.0e-12, % c-axis must be along a principle axis
   if abs(ke'*c) < 1.0e-12, % ray perp to c
      ee = c;
   else,
      disp('ke is parallel to c, not sure what to do here');
      keyboard;
   end
end
eo = unitize(cross(c, ko));
if eo'*eo < 1.0e-12, eo = [0 1 0; 0 0 1; 1 0 0]*c; end

hi = unitize(cross(k1,ei));
hs = unitize(cross(kr,es));
hp = unitize(cross(kr,ep));
ho = unitize(cross(ko,eo));
he = unitize(cross(ke,ee));

A = [
   -es'*ts 0 eo'*ts ee'*ts;
   0 -ep'*tp eo'*tp ee'*tp;
   0 -n1*hp'*ts no*ho'*ts ne*he'*ts;
   -n1*hs'*tp 0 no*ho'*tp ne*he'*tp;
];

rhs = [ei'*ts ei'*tp n1*hi'*ts n1*hi'*tp]';

coef = A \ rhs;
rs = coef(1); rp = coef(2); to = coef(3); te = coef(4);

Itot = ([n1 n1 no ne]'.*(coef.^2))' * abs(([kr kr ko se]'*n));
fprintf('check power conservation: %.8e\n',Itot./(k1'*n));
      

return

function a = unitize(c)
b = sqrt(c'*c);
if abs(b)<1.0e-12, a = zeros(size(c));
else a = c./b;
end
return
