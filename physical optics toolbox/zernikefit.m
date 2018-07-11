function [Z, fRes, rmsresidual, R] = zernikefit(varargin)
% [Z, residual, rmsresidual, R] = zernikefit(x,y,f,n,R,poly_order)
%
% x, y, f are matrices of identical size
% or if x,y are vectors, f is matrix size(length(y),length(x))
% f is the data to be fit, f = f(x,y)
% n (optional, default = 36), fit f to first n zernike polynomials, or if n
%    is a vector, it is a list of the Zernike polynomial designations to
%    use.
% R (optional, default = max(abs([x; y])), radius to normalize
%    x,y grid. Data where sqrt(x.^2 + y.^2) > R is not used in the fit.
% poly_order (optional, default = 'ISO','Korechoff'), 'Noll','Zemax' => Noll ordering 
%
% the polynomials are as
% defined in the memo, R.P. Korechoff, "Extended list of Zernike
% Polynomials," Oct. 13, 2004. sim-lib document #22936
% or Noll ordering if specified.
%
% Z = vector of coefficients, Z(i) is for the i-th polynomial
% rmsresidual = rms(A*Z - f)
% residual = A*Z - f
% R = normalization radius used to determine Z

%persistent P;
 
[x, y, f, n, R, fRes, nuse, poly_order] = ValidateInputs(varargin{:});

%if isempty(P),
    P = zernikepolynomials(poly_order);
%end

r = sqrt(x.^2 + y.^2);
t = atan2(y,x);

M = length(x(:));

A = zeros(M,length(n));
rp = r(:)./R;
for ii = 1:length(n),
    A(:,ii) = P{n(ii)}(rp,t(:));
end

Z = A \ f(:);

% residuals
%r = A*Z - f(:);
r = f(:) - A*Z;

rmsresidual = rms(r(:));
normr = norm(r);
fRes(nuse) = r;

clear A;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y,f,n,R,fRes,nuse,poly_order] = ValidateInputs(varargin)
% [Z, residual, rmsresidual, R] = zernikefit(x,y,f,n,R,poly_order)
% validate inputs

MAXN = 36;

if length(varargin) < 3, error('usage: [Z, normr] = zernikefit(x,y,f,n,R)'); end
xp = varargin{1};
yp = varargin{2};
fp = varargin{3};
% either x,y,f are same size, or f = length(x) x length(y)
if ~(all(size(xp)==size(fp)) && all(size(yp)==size(fp))),
    if ~all(size(fp) == [length(yp),length(xp)]),
        error('input dimension mismatch');
    end
    % x and y are vectors, and f is matrix
    [xp, yp] = meshgrid(xp,yp);
else,
    % x, y, f are same size, use what we got
    % eliminate any NaN
    iuse = ~isnan(fp);
    if ~all(iuse),
        xp = xp(iuse);
        yp = yp(iuse);
        fp = fp(iuse);
    end
end

    function n = CheckN
        a = varargin{4};
        if ~isempty(a),
            if isscalar(a),
                if a > MAXN,
                    warning(['requested order too large, using default ' num2str(MAXN)]);
                    n = [1:MAXN];
                elseif a <= 0,
                    error(['invalid order ' num2str(a)]);
                else,
                    n = [1:a];
                end
            else,
                % assume input is a list, do nothing
                n = a;
            end
        else,
            n = [1:MAXN]; % default value
        end
    end

% radius of each point
rp = sqrt(xp(:).^2 + yp(:).^2);

% defaults for the remaining parameters
n = [1:MAXN];
R = max(rp);
nuse = rp >= 0; % in other words, all the points
poly_order = 'polar';

if length(varargin) >= 4, n = CheckN; end
if length(varargin) >= 5,
    R = varargin{5};
    if ~isempty(R),
        nuse = rp <= R;
    else
        R = max(rp);
    end
end
if length(varargin) >= 6, poly_order = varargin{6}; end

x = xp(nuse);
y = yp(nuse);
f = fp(nuse);

% create residual map for output, it should have same form as fp at input
fRes = nan(size(fp));

end