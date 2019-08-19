function A2d = dh2d(A, nx, whmd, varargin)
% A2d = dh2d(A, nx, whmd)
%
% A is 1-d
% [nx, nx] is size of output A2d
% whmd is 1-d, 0-offset, length(whmd) = length(A)
%
% see hcim/efc/propagate.py:
%    self.whmd
%    def dh2d

if length(A) == 2*length(whmd)
    % assume unraveled complex
    A = A(1:2:end) + 1j*A(2:2:end);
end

A2d = nan(nx);
A2d(whmd+1) = A; % 1-offset, whereas python is 0-offset
A2d = fftshift(A2d);
A2d = A2d.';     % python is row-major

