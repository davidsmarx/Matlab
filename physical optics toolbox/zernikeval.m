function f = zernikeval(Z,x,y,R,poly_order)
% f = zernikeval(Z,x,y,R)
% f = zernikeval(Z,x,y,R,poly_order)
%
% Z is a vector of zernike coefficients
% Z[1] = Z0, Z[2] = Z1, ... as returned by zernikefit()
% x, y are matrices of identical size
% R (optional, default = max(sqrt(x^2+y^2)), normalization radius.
% The normalization radius must be the same that was used to generate the
% coefficients, Z, in zernikefit().
%
% poly_order (optional, default = 'ISO',
%        defined in the memo, R.P. Korechoff, "Extended list of Zernike
%        Polynomials," Oct. 13, 2004. sim-lib document #22936)
%    'Noll', or 'Zemax' uses Noll ordering
%    'Laplacian' uses the Laplacian of the zernike polynomials
%

if exist('poly_order','var'),
    switch lower(poly_order),
        case 'laplacian'
            P = zernikepolynomialsLaplacian;
        case {'iso'},
            P = zernikepolynomials('polar');
        case {'noll','zemax'},
            P = zernikepolynomials('noll');
        case 'codevfringe'
            P = zernikepolynomials('codevfringe');
        otherwise
            error(['unknown polynomial ordering: ' poly_order]);
    end
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

