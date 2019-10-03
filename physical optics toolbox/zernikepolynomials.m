function P = zernikepolynomials(coords)
% P = zernikepolynomials(coords)
%
% coords (optional, default = 'polar') = 'rect'angular or 'polar' or 'noll'
% for Noll (Zemax) ordering
%
% P = cell array of anonymous functions, each a normalized polynomial
% function of (x,y) or (r,t). For example, P{4}(r,t) = sqrt(6)*r^2*sin(2*t)
% is the fourth zernike polynomial. The polynomials are as listed in,
% R.P. Korechoff, "Extended list of Zernike polynomials," JPL Interoffice
% Memorandum, Oct. 13, 2004.
%
% Publications which follow ISO standard 010110
% Cubalchini (1979)
% Malacara 1978, 1992
% Mahajan 1994
% Kaeri & Shannon 1987
% Tyson 1991, 1998
% Love 1997 ****
% Lian, Grimm, Goelz, Bille 1994
% Schwiegerling, Greivenkamp, & Miller 1995
% Liang & Williams 1997

if nargin == 0, coords = 'polar'; end

switch lower(coords),
    case {'pol','polar','p'},
        P = {
            @(r,t) 1                 % zemax term 1, p. 209 of manual
            @(r,t) 2.*r.*sin(t)      % zemax term 3
            @(r,t) 2.*r.*cos(t)      % zemax term 2
            @(r,t) sqrt(6).*r.^2.*sin(2*t)  % zemax term 5
            @(r,t) sqrt(3).*(2.*r.^2 - 1) % 5 zemax term 4
            @(r,t) sqrt(6).*r.^2.*cos(2*t) % 6 zemax term 6
            @(r,t) sqrt(8).*r.^3.*sin(3*t) % 7 zemax term 9
            @(r,t) sqrt(8).*(3*r.^3 - 2*r).*sin(t) % 8, zemax term 7
            @(r,t) sqrt(8).*(3*r.^3 - 2*r).*cos(t) % 9, zemax term 8
            @(r,t) sqrt(8).*r.^3.*cos(3*t)  % 10, zemax term 10
            @(r,t) sqrt(10).*r.^4.*sin(4*t) % 11, zemax term 15
            @(r,t) sqrt(10).*(4*r.^4 - 3*r.^2).*sin(2*t) % 12, zemax term 13
            @(r,t) sqrt(5).*(6*r.^4 - 6*r.^2 + 1) % 13, zemax term 11
            @(r,t) sqrt(10).*(4*r.^4 - 3*r.^2).*cos(2*t) % 14, zemax term 12
            @(r,t) sqrt(10).*r.^4.*cos(4*t) % 15, zemax term 14
            @(r,t) 2*sqrt(3).*r.^5.*sin(5*t)
            @(r,t) 2*sqrt(3).*(5*r.^5 - 4*r.^3).*sin(3*t)
            @(r,t) 2*sqrt(3).*(10*r.^5 - 12*r.^3 + 3*r).*sin(t)
            @(r,t) 2*sqrt(3).*(10*r.^5 - 12*r.^3 + 3*r).*cos(t)
            @(r,t) 2*sqrt(3).*(5*r.^5 - 4*r.^3).*cos(3*t)
            @(r,t) 2*sqrt(3)*r.^5.*cos(5*t) % 21
            @(r,t) sqrt(14)*r.^6.*sin(6*t)
            @(r,t) sqrt(14)*(6*r.^6 - 5*r.^4).*sin(4*t)
            @(r,t) sqrt(14)*(15*r.^6 - 20*r.^4 + 6*r.^2).*sin(2*t)
            @(r,t) sqrt(7)*(20*r.^6 - 30*r.^4 + 12*r.^2 - 1) % 25
            @(r,t) sqrt(14)*(15*r.^6 - 20*r.^4 + 6*r.^2).*cos(2*t)
            @(r,t) sqrt(14)*(6*r.^6 - 5*r.^4).*cos(4*t) % 27
            @(r,t) sqrt(14)*r.^6.*cos(6*t)
            @(r,t) 4*r.^7.*sin(7*t)
            @(r,t) 4*(7*r.^7 - 6*r.^5).*sin(5*t)
            @(r,t) 4*(21*r.^7 - 30*r.^5 + 10*r.^3).*sin(3*t)
            @(r,t) 4*(35*r.^7 - 60*r.^5 + 30*r.^3 - 4*r).*sin(t)
            @(r,t) 4*(35*r.^7 - 60*r.^5 + 30*r.^3 - 4*r).*cos(t)
            @(r,t) 4*(21*r.^7 - 30*r.^5 + 10*r.^3).*cos(3*t)
            @(r,t) 4*(7*r.^7 - 6*r.^5).*cos(5*t)
            @(r,t) 4*r.^7.*cos(7*t) % 36
            };
    case {'noll','zemax'},
        P = {
            @(r,t) 1                 % zemax term 1, p. 209 of manual
            @(r,t) 2.*r.*cos(t)      % zemax term 2
            @(r,t) 2.*r.*sin(t)      % zemax term 3
            @(r,t) sqrt(3).*(2.*r.^2 - 1) % 5 zemax term 4
            @(r,t) sqrt(6).*r.^2.*sin(2*t)  % zemax term 5
            @(r,t) sqrt(6).*r.^2.*cos(2*t) % 6 zemax term 6
            @(r,t) sqrt(8).*(3*r.^3 - 2*r).*sin(t) % 8, zemax term 7
            @(r,t) sqrt(8).*(3*r.^3 - 2*r).*cos(t) % 9, zemax term 8
            @(r,t) sqrt(8).*r.^3.*sin(3*t) % 7 coma in y, zemax term 9
            @(r,t) sqrt(8).*r.^3.*cos(3*t)  % 10, zemax term 10
            @(r,t) sqrt(5).*(6*r.^4 - 6*r.^2 + 1) % 13, zemax term 11
            @(r,t) sqrt(10).*(4*r.^4 - 3*r.^2).*cos(2*t) % 14, zemax term 12
            @(r,t) sqrt(10).*(4*r.^4 - 3*r.^2).*sin(2*t) % 12, zemax term 13
            @(r,t) sqrt(10).*r.^4.*cos(4*t) % 15, zemax term 14
            @(r,t) sqrt(10).*r.^4.*sin(4*t) % 11, zemax term 15
            @(r,t) 2*sqrt(3).*(10*r.^5 - 12*r.^3 + 3*r).*cos(t) % zemax 16
            @(r,t) 2*sqrt(3).*(10*r.^5 - 12*r.^3 + 3*r).*sin(t) % zemax 17
            @(r,t) 2*sqrt(3).*(5*r.^5 - 4*r.^3).*cos(3*t) % zemax 18
            @(r,t) 2*sqrt(3).*(5*r.^5 - 4*r.^3).*sin(3*t) % zemax 19
            @(r,t) 2*sqrt(3)*r.^5.*cos(5*t) % zemax 20
            @(r,t) 2*sqrt(3).*r.^5.*sin(5*t) % zemax 21
            @(r,t) sqrt(7)*(20*r.^6 - 30*r.^4 + 12*r.^2 - 1) % zemax 22
            @(r,t) sqrt(14)*(15*r.^6 - 20*r.^4 + 6*r.^2).*sin(2*t) % zemax 23
            @(r,t) sqrt(14)*(15*r.^6 - 20*r.^4 + 6*r.^2).*cos(2*t) % zemax 24
            @(r,t) sqrt(14)*(6*r.^6 - 5*r.^4).*sin(4*t) % zemax 25
            @(r,t) sqrt(14)*(6*r.^6 - 5*r.^4).*cos(4*t) % zemax 26
            @(r,t) sqrt(14)*r.^6.*sin(6*t) % zemax 27
            @(r,t) sqrt(14)*r.^6.*cos(6*t) % zemax 28
            @(r,t) 4*(35*r.^7 - 60*r.^5 + 30*r.^3 - 4*r).*sin(t) % zemax 29
            @(r,t) 4*(35*r.^7 - 60*r.^5 + 30*r.^3 - 4*r).*cos(t) % zemax 30
            @(r,t) 4*(21*r.^7 - 30*r.^5 + 10*r.^3).*sin(3*t) % zemax 31
            @(r,t) 4*(21*r.^7 - 30*r.^5 + 10*r.^3).*cos(3*t) % zemax 32
            @(r,t) 4*(7*r.^7 - 6*r.^5).*sin(5*t) % zemax 33
            @(r,t) 4*(7*r.^7 - 6*r.^5).*cos(5*t) % zemax 34
            @(r,t) 4*r.^7.*sin(7*t) % zemax 35
            @(r,t) 4*r.^7.*cos(7*t) % zemax 36
            };
        
    case 'codevfringe'
        P = {
           @(r,t) 1 %1 Piston (constant)
           @(r,t) r.*cos(t) % ? Distortion - Tilt (x-axis) 2
           @(r,t) r.*sin(t) % ? Distortion - Tilt (y-axis) 3
           @(r,t) 2.*r.^2 - 1 % Defocus - Field curvature 4
           @(r,t) r.^2.*cos(2*t) % Astigmatism, Primary (axis at 0° or 90°) 5
           @(r,t) r.^2.*sin(2*t) % Astigmatism, Primary (axis at ±45°) 6
           @(r,t) (3*r.^3 - 2*r) .* cos(t) % ? Coma, Primary (x-axis) 7
           @(r,t) (3*r.^3 - 2*r) .* sin(t) %? Coma, Primary (y-axis) 8
           @(r,t) 6*r.^4 - 6*r.^2 + 1 % Spherical Aberration, Primary 9
           @(r,t) r.^3 .* cos(3*t) % Trefoil, Primary (x-axis) 10
           @(r,t) r.^3 .* sin(3*t) % Trefoil, Primary (y-axis) 11
           @(r,t) (4*r.^4 - 3*r.^2) .* cos(2*t) % Astigmatism, Secondary (axis at 0° or 90°) 12
           @(r,t) (4*r.^4 - 3*r.^2) .* sin(2*t) % Astigmatism, Secondary (axis at ±45°) 13
           @(r,t) (10*r.^5 - 12*r.^3 + 3*r) .* cos(t) % Coma, Secondary (x-axis) 14
           @(r,t) (10*r.^5 - 12*r.^3 + 3*r) .* sin(t) % Coma, Secondary (y-axis) 15
           @(r,t) 20*r.^6 - 30*r.^4 + 12*r.^2 - 1 % Spherical Aberration, Secondary 16
           @(r,t) r.^4 .* cos(4*t) % Tetrafoil, Primary (x-axis) 17
           @(r,t) r.^4 .* sin(4*t) % Tetrafoil, Primary (y-axis) 18
           @(r,t) (5*r.^5 - 4*r.^3) .* cos(3*t) % Trefoil, Secondary (x-axis) 19
           @(r,t) (5*r.^5 - 4*r.^3) .* sin(3*t) % Trefoil, Secondary (y-axis) 20
           @(r,t) (15*r.^6 - 20*r.^4 + 6*r.^2) .* cos(2*t) % Astigmatism, Tertiary (axis at 0° or 90°) 21
           @(r,t) (15*r.^6 - 20*r.^4 + 6*r.^2) .* sin(2*t) % Astigmatism, Tertiary (axis at ±45°) 22
           @(r,t) (35*r.^7 - 60*r.^5 + 30*r.^3 - 4*r) .* cos(t) % Coma, Tertiary (x-axis) 23
           @(r,t) (35*r.^7 - 60*r.^5 + 30*r.^3 - 4*r) .* sin(t) % Coma, Tertiar (y-axis) 24
           @(r,t) 70*r.^8 - 140*r.^6 + 90*r.^4 - 20*r.^2 + 1 % Spherical Aberration, Tertiary 25
           @(r,t) r.^5 .* cos(5*t) % Pentafoil, Primary (x-axis) 26
           @(r,t) r.^5 .* sin(5*t) % Pentafoil, Primary (y-axis) 27
           @(r,t) (6*r.^6 - 5*r.^4) .* cos(4*t) % Tetrafoil, Secondary (x-axis) 28
           @(r,t) (6*r.^6 - 5*r.^4) .* sin(4*t) % Tetrafoil, Secondary (y-axis) 29
           @(r,t) (21*r.^7 - 30*r.^5 + 10*r.^3) .* cos(3*t) % Trefoil, Tertiary (x-axis) 30
           @(r,t) (21*r.^7 - 30*r.^5 + 10*r.^3) .* sin(3*t) % Trefoil, Tertiary (y-axis) 31
           @(r,t) (56*r.^8 - 105*r.^6 + 60*r.^4 - 10*r.^2) .* cos(2*t) % Astigmatism, Quaternary (axis at 0° or 90°) 32
           @(r,t) (56*r.^8 - 105*r.^6 + 60*r.^4 - 10*r.^2) .* sin(2*t) % Astigmatism, Quaternary (axis at ±45°) 33
           @(r,t) (126*r.^9 - 280*r.^7 + 210*r.^5 - 60*r.^3 + 5*r) .* cos(t) % Coma, Quaternary (x-axis) 34
           @(r,t) (126*r.^9 - 280*r.^7 + 210*r.^5 - 60*r.^3 + 5*r) .* sin(t) % Coma, Quaternary (y-axis) 35
           @(r,t) 252*r.^10- 630*r.^8 + 560*r.^6 - 210*r.^4 + 30*r.^2 - 1 % Spherical Aberration, Quaternary 36
           @(r,t) 924*r.^12 - 2772*r.^10 + 3150*r.^8 - 1680*r.^6 + 420*r.^4 - 42*r.^2 + 1 %Spherical Aberration, Quinternary 37
           };

    case {'rect','rectangular','r'},
        P = {
            @(x,y) 1
            @(x,y) 2*y
            @(x,y) 2*x
            @(x,y) sqrt(6).*2.*x.*y
            @(x,y) sqrt(3)*(2*(x.^2 + y.^2) - 1)
            @(x,y) sqrt(6)*(x.^2 - y.^2)
            @(x,y) sqrt(8)*(3*x.^2 - y.^2).*y
            @(x,y) sqrt(8)*(3*(x.^2 + y.^2) - 2).*y
            @(x,y) sqrt(8)*(3*(x.^2 + y.^2) - 2).*x
            @(x,y) sqrt(8)*(x.^2 - 3*y.^2).*x
            @(x,y) sqrt(10)*4*(x.^3.*y - x.*y.^3) % 11, zemax term 15
            @(x,y) sqrt(10)*2*x.*y.*(4*(x.^2 + y.^2) - 3) % 12, zemax term 13
            @(x,y) sqrt(5)*(6* (x.^2 + y.^2).^2 - 6*(x.^2 + y.^2) + 1)  % 13, zemax term 11 : 6*x^4 + 12*x^2*y^2 - 6*x^2 + 6*y^4 - 6*y^2 + 1
            @(x,y) sqrt(10)*(4*(x.^2 + y.^2) - 3).*(x.^2 - y.^2) % 14, zemax term 12 : 4*x^4 - 3*x^2 - 4*y^4 + 3*y^2
            @(x,y) sqrt(10)*(x.^4 - 6.*x.^2.*y.^2 + y.^4) % 15, zemax term 14
            };
    otherwise,
        error('invalid input');
end

        
