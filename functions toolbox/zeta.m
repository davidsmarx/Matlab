function [f] = zeta(z)
%ZETA  Riemann Zeta function
%
%usage: f = zeta(z)
%
%tested on version 5.3.1
%
%      This program calculates the Riemann Zeta function
%      for the elements of Z using the Dirichlet deta function.
%      Z may be complex and any size.
%
%      Has a pole at z=1, zeros for z=(-even integers),
%      infinite number of zeros for z=1/2+i*y
%
%
%see also: Eta, Etan, Lambda, Betad, Bern, Euler
%see also: mhelp zeta

%Paul Godfrey
%pgodfrey@intersil.com
%3-24-01

zz=2.^z;
k = zz./(zz-2);

f=k.*deta(z,1);

p=find(z==1);
if ~isempty(p)
   f(p)=Inf;
end

return
