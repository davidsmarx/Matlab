function MTF=MTF(Lambda, F, f)
% MTF(Wavelengths, Focal_Ratio, Frequencies)
%
% Returns the diffraction limited monochromatic MTF for each of the wavelengths given (rows)
% and for each of the frequencies given (columns).
% 
% If the frequencies are given in cycles per millimetre, the wavelengths should also be in mm.
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


% Copyright 2003, Defencetek, CSIR
% $Revision: 1.1 $

% First ensure wavelengths are a column vector
els = size(Lambda,1) * size(Lambda,2);
if (els == size(Lambda,2))
    Lambda = Lambda';
end
if (els ~= size(Lambda,1))
    Lambda = Lambda(:,1);
end

% Next ensure frequencies are a row vector
els = size(f,1) * size(f,2);
if (els == size(f,1))
    f=f';
end
if (els ~= size(f,2))
    f=f(:,1);
end

% Ensure F scalar, positive
F = abs(F(1));

% Now mesh the wavelengths and the frequencies
[x,y]=meshgrid(f, Lambda);
% Find frequencies above the cutoff

cutoffs = ones(size(Lambda,1),1)./(F * Lambda);
% For each row, set any frequencies exceeding cutoff to cutoff
for i=1:(size(Lambda,1))
    j = find(x(i,:) > cutoffs(i));
    x(i,j) = cutoffs(i);
end


% Compute Lambda * frequency * Focal_ratio for each site in the matrix
phi = acos(F * x.*y);
csphi = cos(phi) .* sin(phi);
MTF = 2.0 * (phi - csphi) / pi;


