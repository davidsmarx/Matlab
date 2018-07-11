function U1 = ApplyLens(U0,rn,lens,lam)
% U1 = ApplyLens(U0,rn,lens,lam)
%
% lens.f = focal length
% lens.D = clear aperture diameter

U1 = zeros(size(U0));

n0 = rn < (lens.D/2);

phi = (pi./lam./lens.f).*(rn(n0).^2);
U1(n0) = U0(n0).*exp(-j*phi);

