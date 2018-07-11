function f = zernikeval(Z,x,y,R,uselaplacian)
% f = zernikeval(Z,x,y,R)
% f = zernikeval(Z,x,y,R,uselaplacian)
%
% Z is a vector of zernike coefficients
% Z[1] = Z0, Z[2] = Z1, ... as returned by zernikefit()
% x, y are matrices of identical size
% R (optional, default = max(sqrt(x^2+y^2)), normalization radius.
% The normalization radius must be the same that was used to generate the
% coefficients, Z, in zernikefit().
%
% The polynomials are as
% defined in the memo, R.P. Korechoff, "Extended list of Zernike
% Polynomials," Oct. 13, 2004. sim-lib document #22936
%
% Z = vector of coefficients, Z(i) is for the i-th polynomial
%
% if uselaplacian == true (default = false) then the Laplacian of the
% Zernike polynomials are used.

if exist('uselaplacian','var') & uselaplacian == true,
    P = zernikepolynomialsLaplacian;
else,
    P = zernikepolynomials('polar');
end

if ~exist('R','var'),
    R = max(r(:)) % max radius to normalize coordinates
end

r = sqrt(x.^2 + y.^2);
t = atan2(y,x);

N = length(Z);    % first N polynomials are used

f = zeros(size(r));
rp = r./R;
for ii = 1:N,
    f = f + Z(ii).*P{ii}(rp,t);
end
