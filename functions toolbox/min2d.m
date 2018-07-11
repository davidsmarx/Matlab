function [amax, row, col] = min2d(A)
% [min row col] = min2d(A)
% input A is a 2-dimensional matrix
% min is the minimum of abs(A)
% row, col is the location of the maximum

if nargin == 0, error('usage: [amax, row, col] = min2d(A)'); end

if ~isreal(A), [bmax imax] = min(abs(A));
else [bmax imax] = min(A);
end

[bmaxmax imaxmax] = min(bmax);

amax = bmaxmax;
row = imax(imaxmax);
col = imaxmax;

return

