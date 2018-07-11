function [r, Ar] = RadialMean(varargin)
% [r, Ar] = RadialMean(A)
% [r, Ar] = RadialMean(x, A) % assumes y == x
% [r, Ar] = RadialMean(x, y, A)
% [r, Ar] = RadialMean(..., N) % N = length of return r (default = 128)

% D. Marx 2018-05-18 initial version

% default
Nr = 128;

IsMatrix = @(a) ~isscalar(a) && ~isvector(a);
IsVector = @(a) ~isscalar(a) && min(size(a)) == 1;

switch nargin,
    case 1,
        A = varargin{1};
        [x, y, X, Y, R] = CreateGrid(A);
        
    case 2,
        if IsMatrix(varargin{1}),
            [A, Nr] = deal(varargin{:});
            [x, y, X, Y, R] = CreateGrid(A);
        elseif IsVector(varargin{1})
            [x, A] = deal(varargin{:});            
            y = x;
            [X, Y] = meshgrid(x, y); R = sqrt(X.^2 + Y.^2);
        else % isscalar(varargin{1})
            [dx, A] = deal(varargin{:});
            [x, y, X, Y, R] = CreateGrid(A, dx);
        end
        
    case 3,
        if IsMatrix(varargin{2}),
            [x, A, Nr] = deal(varargin{:});
            y = x;
        else
            [x, y, A] = deal(varargin{:});
        end
        [X, Y] = meshgrid(x, y); R = sqrt(X.^2 + Y.^2);
        
    case 4,
        [x, y, A, Nr] = deal(varargin{:});
        [X, Y] = meshgrid(x, y); R = sqrt(X.^2 + Y.^2);
        
end

rmax = max(abs([x(:); y(:)]));
dr   = rmax./Nr;

r = linspace(dr/2,rmax-dr/2,Nr).';
Ar = zeros(size(r));

for ir = 1:Nr,
    Ar(ir) = mean( A(R>r(ir)-dr/2 & R<=r(ir)+dr/2) );
end

