function f = zernikeval(Z,x,y,R,varargin)
% f = zernikeval(Z,x,y,R)
% f = zernikeval(Z,x,y,R,poly_order)
% f = zernikeval(Z,x,y,R,..,options)
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
% options:
%    nz: list of Zernike #'s corresponding to each coeff in Z, same size as Z
%        default: nz = 1:length(Z)

NZ = CheckOption('nz', 1:length(Z), varargin{:});
if ~isequal(length(Z), length(NZ)),
    error('list of Z not same size as list of NZ');
end
NZ = reshape(NZ, size(Z));

if ~isempty(varargin),
    poly_order = varargin{1};
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
            fprintf('using default poly order: Noll\n');
            P = zernikepolynomials('noll');
    end
else,
    P = zernikepolynomials('polar');
end

if ~exist('R','var'),
    R = max(r(:)) % max radius to normalize coordinates
end

if ~isequal(length(Z), length(NZ)),
    error('list of Z not same size as list of NZ');
end

r = sqrt(x.^2 + y.^2);
t = atan2(y,x);

f = zeros(size(r));
rp = r./R;
for ii = 1:length(Z)
    pp = P{NZ(ii)};
    f = f + Z(ii).*pp(rp,t);
end

