function [sp,se,st,su] = pdesd(p,e,t,u,sdl)
% [sp,se,st,su] = psesd(p,e,t,u,sdl)
% return the p,e,t and solution u values within the subdomain sdl

[ip, cp] = pdesdp(p,e,t,sdl);
sp = [p(:,ip), p(:,cp)];

it = pdesdt(t,sdl); st = t(:,it);

ie = pdesde(e,sdl); se = e(:,ie);

if length(u) == length(p) % then u is point data
   su = [u(ip); u(cp)];
elseif length(u) == length(t) % then u is triangle node data
   su = u(it);
else
   error('data u does not match point or node data')
end
      