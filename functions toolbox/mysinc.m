function y=sinc(x)
% SINC  SINC(X) is SIN(X)/X  with the singularity at zero removed.

%       Drea Thomas     1992

if isempty(x), error('input is empty'); end
z = x ~= 0;
y=x;
y(z) = sin(x(z))./x(z);
%y(~z)=ones(sum(sum(~z)),1);
y(~z) = 1.0;
