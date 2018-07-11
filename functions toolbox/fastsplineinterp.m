function yi = fastsplineinterp(y,xi)
% yi = fastsplineinterp(y,xi)
% fast cubic spline of the vector y at the points xi
% y is a vector of samples taken on a uniformly spaced grid 1..Ny
% xi is a list of fractional grid values between 1 and Ny
% if xi is a scalar between 0 and 1, then it is interpreted as a shift value,
% and yi is y shifted by xi. Otherwise, yi is calculated for each xi.
% yi are the interpolated values of y on the points xi.

z1 = -2 + sqrt(3);

if nargin ~= 2, error('usage: yi = fastsplineinterp(y,xi)'); end

y  = y(:); % make sure y is a column vector
Nx = length(xi);
Ny = length(y);

ap = [1 -z1]; bp = 1;
cp = filter(bp,ap,y); cp(1) = sum(y.* (z1.^[0:Ny-1]') );

cpp= cp(Ny:-1:1); % time reverse cp for anti-causal filter
am = [1 -z1]; bm = -z1;
cm = filter(bm,am,cpp); cm(1) = z1*(cp(1) + z1*cp(2))./(1-z1.^2);
cm = cm(Ny:-1:1);

if (Nx == 1) & (xi == 0),
   yi = y;
elseif (Nx == 1) & (xi > 0) & (xi < 1),
   b3 = B3(xi + [-2 -1 0 1 2]);
   yi = 6*filter(b3,1,cm);
   % but yi is now delayed by one sample
   yi = [yi(3:Ny); 0; 0];
elseif Nx==Ny,
   yi = y;
   k = [-2:2]';
   for n = 3:Ny-2,
      yi(n) = 6*sum(cm(n+k).*B3((xi(n)-n)-k));
   end
else,
   warning('this option not implemented');
end



return


function y = B3(x)

ax = abs(x);
y = zeros(size(ax));
x1 = ax<1;
y(x1) = (2/3) - ax(x1).^2 + 0.5*(ax(x1).^3);
x2 = ( (ax<2)&(ax>=1) );
y(x2) = ( (2 - ax(x2)).^3 )/6;

return
