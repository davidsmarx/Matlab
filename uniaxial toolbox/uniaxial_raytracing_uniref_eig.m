function [ke, se] = uniaxial_raytracing_uniref_eig(k1, n, c, mo, me)
% from Beyerle & McDermid paper
% k1, n, c are unit vectors of incident phase velocity, surface normal, axis orientation
% m1 = index of incident medium
% mo, me = ord. & ext. index of crystal
% ke, se = extra-ordinary phase and ray unit vector directions
% m2 = index of refraction for extraordinary
% ko = ordinary phase and ray direction

% make sure vectors are unit
k1 = k1(:)./mag(k1);
n  = n(:) ./mag(n);
c  = c(:) ./mag(c);

% form gamma and orthogonal transform R
phi = atan2(c(2),c(1)); 
cph = cos(phi); sph = sin(phi);
cth = c(3);     sth = sqrt(1-cth.^2);
R = [cth*cph cth*sph -sth; -sph cph 0; sth*cph sth*sph cth];
einv = R' * diag([1./mo^2 1./mo^2 1./me^2]) * R;
g = R' * diag([1./me 1./me 1./mo]) * R;
% determine m1 for k1
m1 = 1./mag(g*k1);

% initialize
m2 = me; mold = mo;
while abs(m2-mold) > 1.0e-6,
   mold = m2;
   ke = (m1/m2)*( k1 + (-k1'*n - sqrt( (m2/m1)^2 - 1 + (k1'*n)^2 ))*n );
   % form matrix k x k x
   K = cross([ke ke ke],eye(3));
   [V, D] = eig(einv*(K*K));
   % choose the appropriate e.v. e.vectors, V, are in direction of E-field
   [atmp, iv] = max(abs(V'*(K*K*c)));
   m2 = 1/sqrt(-D(iv,iv));
   keyboard;
end

se = einv \ ke;  % from Eq. 7 of paper
se = se./mag(se);

% check fresnel's equation:
kk = R*ke; % move ke into crystal coordinates
checkfresn = abs( sum((kk.^2)./([me me mo]'.^2)) - 1./(m2^2) );
if checkfresn > 1.0e-9, 
   warning(['warning: extra-ordinary does not satisfy fresnel condition: '...
         num2str(checkfresn,'%.2e')]);
end

% ordinary direction (Snell's Law)
sqterm = sqrt( (mo/m1)^2 -1 + (k1'*n)^2 );
ko = (m1/mo)*(k1 + n*(sqterm - k1'*n));
ko = ko/sqrt(ko'*ko);

return

function a = mag(c)
a = sqrt(c'*c);
return
 


