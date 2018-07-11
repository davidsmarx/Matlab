function n_air = n_air(lambda, T, P)
% n_air(lambda, T, P) returns the refractive index of air computed using the same formula used by ZEMAX
% See the section on Index of Refraction Computation in the Thermal Analysis chapter of the ZEMAX manual.
% 
% lambda is wavelength in microns.
% T is temperature in Celsius.
% P is relative air pressure.
% This function returns a matrix with lambda varying from row to row, temperature varying from column to column
% and pressure varying in the depth dimension.
%

% MZDDE - The ZEMAX DDE Toolbox for Matlab.
% Copyright (C) 2002-2004 Defencetek, CSIR
% Contact : dgriffith@csir.co.za
% 
% This file is part of MZDDE.
% 
%  MZDDE is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  MZDDE is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with MZDDE (COPYING.html); if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%


% Copyright 2002, Defencetek, CSIR
% $Revision: 1.2 $

% Meshgrid the data 
[T, lambda, P] = meshgrid(T, lambda, P);

% Compute the reference refractive indices across the wavelengths
n_ref = 1 + (6432.8 + (2949810 * lambda.^2) ./ (146 * lambda.^2 - 1) + (25540 * lambda.^2) ./ (41 * lambda.^2 - 1)) * 1e-8;

% Finally, compute the full (potentially 3D) data set
n_air = 1 + ((n_ref - 1) .* P) ./ (1 + (T - 15) * 3.4785e-3);
