function [r, Ar] = RadialMean(x,y,A,varargin)
% [r, Ar] = RadialMean(x, y, A, varargin)
% 
% x,  y are vectors defining a 2-d regular grid
% A is size length(y) x length(x)
% 
% options:
%    CheckOption('Nr', 128
%    CheckOption('bMask', true(size(A))
%    CheckOption('thetastartstop', [], options{:});
%
% D. Marx 2018-05-18 initial version
%
% 2021-11-24
%   substantial changes to include bMask

% default
Nr = CheckOption('Nr', 128, varargin{:});
bMask = CheckOption('bMask', true(size(A)), varargin{:});
thetastartstop = CheckOption('thetastartstop', [], varargin{:});

[X, Y] = meshgrid(x,y);
R = hypot(X,Y);
T = atan2(Y,X);

rmax = max(abs([x(:); y(:)]));
dr   = rmax./Nr;

r = linspace(dr/2,rmax-dr/2,Nr).';
Ar = zeros(size(r));

for ir = 1:Nr,
    iuse = R>r(ir)-dr/2 & R<=r(ir)+dr/2 & bMask;
    if ~any(iuse(:)), continue, end

    if ~isempty(thetastartstop),
        iuse = iuse & mod2pi(T - thetastartstop(1)) >= 0 & mod2pi(T - thetastartstop(2)) <= 0;
        %figure, imageschcit(x, y, iuse)
    end
    Ar(ir) = mean( A(iuse) );
end

