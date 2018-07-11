function WriteZemaxWaveFile(File, WaveMatrix, Primary)
% WriteZemaxWaveFile - Writes a list of wavelengths formatted for ZEMAX
%
% Usage : WriteZemaxWaveFile(File, WaveMatrix, Primary)
% File is the file to write.
% WaveMatrix is a matrix of wavelengths and weights (one each per row).
% Primary is the primary wavelength number.
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


% Copyright 2004, Defencetek, CSIR
% $Revision: 1.1 $


fid = fopen(File, 'w');
if (fid == -1)
    disp('Cannot Open File');
    return;
end

fprintf(fid, '%i\r\n', Primary);
for i = 1:size(WaveMatrix,1)
    fprintf(fid, '%f %f\r\n', WaveMatrix(i,1), WaveMatrix(i,2));
end
fclose(fid);
