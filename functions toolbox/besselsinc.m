function y=besselsinc(x)
% besselj(1,x)/x  with the singularity at zero removed.

%       Drea Thomas     1992

if isempty(x), error('input is empty'); end
z = x ~= 0;
y=x;
y(z) = besselj(1,x(z))./x(z);
%y(~z)=ones(sum(sum(~z)),1);
y(~z) = 0.5;
