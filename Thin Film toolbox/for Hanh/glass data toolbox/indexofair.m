function [n, dndt] = indexofair(ll,tt,pp)
% [n, dndt] = indexofair(wavelength,temperature,pressure)
% wavelength = [m]
% temperature = [deg C]
% pressure = [relative pressure, atmosphere, 1.0 = standard atmospheric pressure]

if nargin == 0, error('usage: n = indexofair(wavelength,temperature,pressure)'); end

ll = ll./CConstants.UM; % equations below are for wavelength in microns

syms lam p t real;

l2 = lam^2;
nref = 1 + 1e-8*(6432.8 + 2949810*l2/(146*l2-1) + 25540*l2/(41*l2-1));
nair = 1 + (nref - 1)*p/(1 + 3.4785e-3*(t-15));

[L, T, P] = meshgrid(ll,tt,pp);

n  = squeeze(double(subs(nair,{lam,p,t},{L,P,T})));

if nargout > 1,
   dnairdt = diff(nair,t);
   dndt = squeeze(double(subs(dnairdt,{lam,p,t},{L,P,T})));
end

return

fprintf('index of air = %f\n',n);

dndt = diff(nair,T);
dndp = diff(nair,P);
n1 = double(subs(dndt,{lam,P,T},{ll,pp,tt}))
n2 = double(subs(dndp,{lam,P,T},{ll,pp,tt}))
