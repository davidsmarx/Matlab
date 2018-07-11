function y = Fr_exp(x)
% y = Fr_exp(x)
% the fresnel integral: y = int(exp(j(pi/2) t^2),t,0,x)
% calculated using the rational approximations from
% Abramowitz and Stegun, p.302, 7.3

f = (1 + 0.926*x)./(2+x.*(1.792+x.*3.104));
g = 1./(2+x.*(4.142+x.*(3.492+x.*6.670)));

y = 0.5*(1+j) - exp(j*0.5*pi*x.*x) .* (g + j*f);
