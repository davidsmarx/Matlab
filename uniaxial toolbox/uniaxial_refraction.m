function [eo, ee, es, ep, to, te, rs, rp] = uniaxial_refraction(k1,ko,ke,se,n,m,c,ei)
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
ts = cross(k1,n); if ts'*ts < 1.0e-12, ts = [0 1 0;-1 0 0;0 0 1]*n; else ts = ts./sqrt(ts'*ts); end
tp = cross(n,ts); if abs( (tp'*tp)-1.0 ) > 1.0e-12, error('something wrong with tp'); end

es = ts;
es = cross(kr,n); if es'*es < 1.0e-12, es = -ts; else es = es./sqrt(es'*es); end
ep = unitize(cross(kr,es));
ee = unitize(cross(se, cross(se,ke)));
eo = unitize(cross(c, ko));

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

fprintf('check power conservation: %.8e\n',([n1 n1 no ne]'.*coef)'*coef);

      

return

er = ei - 2*(n'*ei)*n; % direction of reflected E field
de = unitize(cross(ke, cross(se,ke)));


% H field directions from Maxwell's Eqn H = curl(E) = k x E
hi = cross(k1,ei);
hr = cross(kr,er);

% calculate transmission coefficients
A = [-er + (n'*er)*n, eo - (n'*eo)*n, ee - (n'*ee)*n];
rhse = ei - (n'*ei)*n;
B = [-(hr - (n'*hr)*n), (ho - (n'*ho)*n), (he - (n'*he)*n)]*diag(m);
rhsh = m(1)*(hi - (n'*hi)*n);
C = [A; B];
rhs = [rhse; rhsh];

% use svd to solve these 6 linearly dependent eqns to solve for 3 unknowns
[U, S, V] = svds(C);
lam = diag(S); ii = abs(lam)>1.0e-12; lam(ii) = 1./lam(ii); IS = diag(lam);
coef = (V*IS*(U'))*rhs;
%fprintf('check power conservation: %.8e\n',(m.*coef)'*coef);
% now try slash
% coef = C \ rhs;
fprintf('check power conservation: %.8e\n',(m.*coef)'*coef);

% coefficients:
r = coef(1); to = coef(2); te = coef(3);

return

function a = unitize(c)
b = sqrt(c'*c);
if abs(b)<1.0e-12, a = zeros(size(c));
else a = c./b;
end
return
