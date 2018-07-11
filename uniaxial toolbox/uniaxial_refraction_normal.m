function [thee, tht, n] = uniaxial_refraction(thi,phi,n1,no,ne)
% [thee, thet, n] = uniaxial_refraction(thi,phi,n1,no,ne)
% thi = angle of incidence
% phi = angle of c axis pi/2 => normal to surface
% n1  = index of incident medium
% no  = ordinary index
% ne  = extraordinary index
% thee= ray angle (electric field)
% thet= wave angle (phase velocity)
% n   = effective index of refraction (phase velocity)

A   = (cos(phi)./no).^2 + (sin(phi)./ne).^2;
B   = 1./ne.^2 - 1./no.^2;
sthi= sin(thi);

% initialize
s2tht = sin(2*thi);
last  = 10;
cnt   = 0;
while mean(abs(last-s2tht))>1e-9 & cnt < 10,
   last = s2tht;
   n2 = 1./(A + B.*s2tht);
   stht = n1.*sthi./sqrt(n2);
   if stht > 1, error('sin(theta_t) > 1'); end
   s2tht = 2.*stht.*sqrt(1-stht.*stht);
   cnt = cnt + 1;
   %fprintf('s2tht = %f last = %f n = %f\n',s2tht,last,sqrt(n2));
end
fprintf('converge after %d reps\n',cnt);
n = sqrt(n2);
tht = asin(stht);
ctht = sqrt(1 - stht.*stht);

U = [cos(phi) sin(phi); -sin(phi) cos(phi)];
G = diag([1./(no.^2) 1./(ne.^2)]);
GG= U' * G * U;

df = size(stht);
E = GG * [stht(:)'; ctht(:)'];
thee = -atan2(-E(1,:),E(2,:)); % ray direction = E x H = GG[sth cth 0] x [0 0 1] = [Ey -Ex 0]
thee = reshape(thee,df);

return