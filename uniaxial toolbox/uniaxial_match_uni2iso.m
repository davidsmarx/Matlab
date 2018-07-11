function [es, ep, eo, ee, ts, tp, ro, re] = ...
   uniaxial_match_uni2iso(k1,k2,ko,ke,se,n,mi,mo,me,m2,c,ei)
%[es, ep, eo, ee, ts, tp, ro, re] = ...
%  uniaxial_match_uni2iso(k1,k2,ke,ko,se,n,mi,mo,me,m2,c,ei)
%
% k1 direction vector of incident ray
% k2, ko, ke, se direction vectors of transmitted, ord., ext. phase, ext. ray
% n direction vector of surface normal
% mi = index for incident k-vector
% mo, me = index for o and e reflected waves (me is for ke from calculation)
% m2 = index of isotropic transmission medium
% c crystal axis direction vector
% ei direction vector of incident polarization
% outputs:
%  eo, ee, es, ep E-field direction vectors
%  ts, tp transmission coefficients
%  ro, re reflection coefficient
% resulting E-fields are:
%   Eo = sqrt(mo)*ro*eo;
%   Ee = sqrt(me)*re*ee;
%   Et = sqrt(m2)*(ts*es + tp*ep);

% surface tangential vectors:
us = cross(k2,n); if us'*us < 1.0e-12, us = [0 1 0;0 0 1;1 0 0]*n; else us = us./sqrt(us'*us); end
up = cross(n,us); if abs( (up'*up)-1.0 ) > 1.0e-12, error('something wrong with up'); end

ei = unitize(ei);
es = us;
ep = up;
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
ho = unitize(cross(ko,eo));
he = unitize(cross(ke,ee));
hs = unitize(cross(k2,es));
hp = unitize(cross(k2,ep));

A = [
   -eo'*us -ee'*us es'*us ep'*us;
   -eo'*up -ee'*up es'*up ep'*up;
   -mo*ho'*us -me*he'*us m2*hs'*us m2*hp'*us;
   -mo*ho'*up -me*he'*up m2*hs'*up m2*hp'*up;
];

rhs = [ei'*us ei'*up mi*hi'*us mi*hi'*up]';
keyboard;

coef = A \ rhs;
ro = coef(1); re = coef(2); ts = coef(3); tp = coef(4);

Itot = ([mo me m2 m2]'.*(coef.^2))' * abs([ko se k2 k2]'*n);
fprintf('check power conservation: %.8e\n',Itot./(mi*k1'*n));
      

return

function a = unitize(c)
b = sqrt(c'*c);
if abs(b)<1.0e-12, a = zeros(size(c));
else a = c./b;
end
return
