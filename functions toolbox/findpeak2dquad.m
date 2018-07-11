function [xp, yp, Fp] = findpeak2dquad(X, Y, F)
% [xp, yp, Fp] = findpeak2dquad(X, Y, F)
%
% X, Y, F are same size, [Ny, Nx]
%
% quadratic fit (6 terms) to F(X,Y)
% return location and value of peak of the quadratic surface (grad F == 0)

if ~isequal(size(X),size(Y)),
    error('X and Y are not equal size');
end
if ~isequal(size(X),size(F)),
    error('X and F are not equal size');
end

A = [ones(size(X(:))) X(:) Y(:) X(:).^2 X(:).*Y(:) Y(:).^2];

w = A \ F(:);

den = 4*w(4)*w(6) - w(5).^2;
xp = (w(3)*w(5) - 2*w(2)*w(6))./den;
yp = (w(2)*w(5) - 2*w(3)*w(4))./den;
Fp = [1 xp yp xp.^2 xp.*yp yp.^2]*w;
