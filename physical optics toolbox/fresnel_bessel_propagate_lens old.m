function u2 = fresnel_bessel_propagate_lens(r2,a0,a1,f,lam)
% u1 = fresnel_bessel_propagate(r2,a0,a1,f,lam)
% a = RADIUS
% f = focal length
% propagation is to distance of 2f behind lens

pi_lamf = pi./lam./f;

u2tmp = zeros(size(r2));
for ii = 1:length(r2),
    u2tmp(ii) =...
        dblquad(@FresnelBesselDblIntegrand,0,a0,0,a1,1e-17,[],pi_lamf,r2(ii));
end

u2 = -2.*(pi_lamf.^2).*exp(j.*0.5.*pi_lamf.*r2.^2).*u2tmp;

%disp('top level counts:');
%disp(cnts);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = FresnelBesselDblIntegrand(r0,r1,pi_lamf,r2)

y = exp(j.*pi_lamf.*(0.5.*r1.^2 + r0.^2))...
    .*besselj(0,2.*pi_lamf.*r0.*r1)...
    .*besselj(0,pi_lamf.*r1.*r2)...
    .*r0.*r1;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = FresnelBesselLensIntegrand(r1,pi_lamf,a0,r2)

u1 = zeros(size(r1));
for ii = 1:length(r1),
    [u1(ii), cnts(ii)] = quadl(@FresnelBesselIntegrand,0,a0,1e-15,[],pi_lamf,r1(ii));
end
%disp(cnts);
y = u1.*exp(-j.*0.5.*pi_lamf.*r1.^2).*besselj(0,pi_lamf.*r2.*r1).*r1;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = FresnelBesselIntegrand(r0,pi_lamz,r1)
% y = fresnel_bessel_propagate(r0,pi_lamz,r1)
% 
% y = exp(j*pi_lamz.*r0.^2).*besselj(0,2.*pi_lamz.*r1.*r0).*r0;

y = exp(j*pi_lamz.*r0.^2).*besselj(0,2.*pi_lamz.*r1.*r0).*r0;

return