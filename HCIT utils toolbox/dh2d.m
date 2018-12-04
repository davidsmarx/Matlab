function A2d = dh2d(A, nx, whmd)
% A2d = dh2d(A, nx, whmd)
%
% A is 1-d
% [nx, nx] is size of output A2d
% whmd is 1-d, length(whmd) = length(A)
%
% see hcim/efc/propagate.py:
%    self.whmd
%    def dh2d

A2d = zeros(nx);
A2d(whmd+1) = A; % 1-offset, whereas python is 0-offset
A2d = fftshift(A2d);
A2d = A2d.';     % python is row-major

