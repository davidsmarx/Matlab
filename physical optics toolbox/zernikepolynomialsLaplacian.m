function P = zernikepolynomialsLaplacian
% P = zernikepolynomialsLaplacian
%
% P = cell array of anonymous functions, each the Laplacian of the
% corresponding Zernike Polynomial
P = {
    @(r,t) 0
    @(r,t) 0
    @(r,t) 0
    @(r,t) 0
    @(r,t) 8*3.^(1/2)
    @(r,t) 0
    @(r,t) 0
    @(r,t) 48*2.^(1/2).*r.*sin(t)
    @(r,t) 48*2.^(1/2).*r.*cos(t)
    @(r,t) 0
    @(r,t) 0
    @(r,t) 48*10.^(1/2).*sin(2.*t)*r.^2
    @(r,t) 24*5.^(1/2).*(4*r.^2-1)
    @(r,t) 48*10.^(1/2).*cos(2*t).*r^2
    @(r,t) 0
    @(r,t) 0
    @(r,t) 160*3.^(1/2)*r.^3.*sin(3*t)
    @(r,t) 96*3.^(1/2).*r.*sin(t).*(5*r^2-2)
    @(r,t) 96*3.^(1/2).*r.*cos(t).*(5.*r^2-2)
    @(r,t) 160*3.^(1/2).*r.^3.*cos(3*t)
    @(r,t) 0
    @(r,t) 0
    @(r,t) 120*14.^(1/2).*r.^4.*sin(4.*t)
    @(r,t) 240*14.^(1/2).*sin(2*t).*r.^2.*(2*r.^2-1)
    @(r,t) 48*7.^(1/2).*(15*r.^4-10.*r.^2+1)
    @(r,t) 240*14.^(1/2).*cos(2*t).*r.^2.*(2*r.^2-1)
    @(r,t) 120*14.^(1/2).*r.^4.*cos(4*t)
    @(r,t) 0
    @(r,t) 0
    @(r,t) 672*r.^5.*sin(5*t)
    @(r,t) 480*r.^3.*sin(3*t).*(7*r.^2-4)
    @(r,t) 960*r.*sin(t).*(7*r.^4-6*r.^2+1)
    @(r,t) 960*r.*cos(t).*(7*r.^4-6*r.^2+1)
    @(r,t) 480*r.^3*cos(3*t).*(7*r.^2-4)
    @(r,t) 672*r.^5.*cos(5*t)
    @(r,t) 0
};