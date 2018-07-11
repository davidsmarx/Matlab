% GAMATH   Translate and evaluate Mathematica's extended gamma function 
%          to/in MATLAB. Main difference is the definition of the
%          incomplete gamma function. In contrast to MATLAB's
%          implementation, but in agreement with Mathematica's, complex 
%          arguments are allowed.
%
%   usage: Same as Gamma[a(,z0(,z1))] in Mathematica, but vector/matrix
%          arguments (in combination with scalars) are allowed.

% Details:
%
% The simple gamma(x) function is defined alike in Mathematica and MATLAB,
%   complex arguments are allowed in MATLAB through use of GAMMAZ.M 
%   instead of GAMMA.M
%
% MATLAB's incomplete gamma function is Abramowitz&Stegun eq. 6.5.1,
%   however the order in the arguments is switched:
% P(a,x) = 1/gamma(a) int_0^infty exp(-t) t^(a-1) dt  =  GAMMAINC(x,a)
%
% Mathematica's incomplete gamma function is A&S eq. 6.5.3
% Gamma(a,x) = int_x^infty  exp(-t) t^(a-1) dt  (a in Complex)
%            = gamma(a)*(1-P(a,x))
%
% Mathematica's 3-argument version is defined as
% Gamma(a,z0,z1) = Gamma(a,z0)-Gamma(a,z1)
%                = gamma(a)*(-P(a,z0)+P(a,z1))

function f = gamath(a,z0,z1);

if nargin == 1,
  f = gammaz(a);
elseif nargin == 2,
  f = gammaz(a).*(1-gammaiz(z0,a));
elseif nargin == 3,
  f = gammaz(a).*(-gammaiz(z0,a) + gammaiz(z1,a));
else
  error('Wrong number of arguments in GAMATH');
end
