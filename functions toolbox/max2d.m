function [amax, row, col] = max2d(A)
% [max row col] = max2d(A)
% input A is a 2-dimensional matrix
% max is the maximum of abs(A)
% row, col is the location of the maximum

if nargin == 0, disp('usage: [amax, row,col] = max2d(A)'); return, end

% if isreal(A), [bmax imax] = max(A); else, [bmax imax] = max(abs(A)); end
% 
% [bmaxmax imaxmax] = max(bmax);
% 
% amax = bmaxmax;
% row = imax(imaxmax);
% col = imaxmax;

if isreal(A), [amax, imax] = max(A(:)); else, [amax, imax] = max(abs(A(:))); end
[row, col] = ind2sub(size(A),imax(1));

if nargout == 0
    disp([amax, row, col]);
end

return

