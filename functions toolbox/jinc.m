function y=jinc(x)
%JINC BesselJ(1,pi*x)/(pi*x) function.
%   JINC(X) returns a matrix whose elements are the sinc of the elements 
%   of X, i.e.
%        y = BesselJ(1,pi*x)/(pi*x)    if x ~= 0
%          = 1/2                   if x == 0
%   where x is an element of the input matrix and y is the resultant
%   output element.
%
%   % Example of a jinc function for a linearly spaced vector:
%   t = linspace(-5,5);
%   y = jinc(t);
%   plot(t,y);
%   xlabel('Time (sec)');ylabel('Amplitude'); title('Sinc Function')
%
%   See also SQUARE, SIN, COS, CHIRP, DIRIC, GAUSPULS, PULSTRAN, RECTPULS,
%   and TRIPULS.

%   Author(s): T. Krauss, 1-14-93
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.7.4.1 $  $Date: 2004/08/10 02:11:27 $
%   revised 2018-07-16, D. Marx, to use bool indexing

% i=find(x==0);                                                              
% x(i)= 1;      % From LS: don't need this is /0 warning is off                           
% y = besselj(1,pi*x)./(pi*x);                                                     
% y(i) = 1/2;   

ii = ~(x==0);
y = zeros(size(x));
y(~ii) = 0.5;
y(ii) = besselj(1,pi*x(ii))./(pi*x(ii));

return
