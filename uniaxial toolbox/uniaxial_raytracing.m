function [ke, se, m2, ko] = uniaxial_raytracing(k1, n, c, m1, mo, me)
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
g = R' * diag([1./me 1./me 1./mo]) * R;

% equation 14 in paper:
sqterm = sqrt( ( (g*n)'*(g*k1) ).^2 - (mag(g*n).*mag(g*k1)).^2 + (mag(g*n)./m1).^2 );
num = (mag(g*n).^2)*k1 + (sqterm - (g*n)'*(g*k1))*n;
ke = num./mag(num);

se = (g*g)*ke; se = se./mag(se);
m2 = 1./mag(g*ke);

% check fresnel's equation:
kk = R*ke; % move ke into crystal coordinates
checkfresn = abs( sum((kk.^2)./([me me mo]'.^2)) - 1./(m2^2) );
if checkfresn > 1.0e-9, 
   warning(['warning: extra-ordinary does not satisfy fresnel condition: '...
         num2str(checkfresn,'%.2e')]);
end

% ordinary direction
sqterm = sqrt( (mo/m1)^2 -1 + (k1'*n)^2 );
ko = (m1/mo)*(k1 + n*(sqterm - k1'*n));
ko = ko/sqrt(ko'*ko);

return

function a = mag(c)
a = sqrt(c'*c);
return
 


