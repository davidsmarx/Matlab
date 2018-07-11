function [w,n,comment] = read_nk(filename,scale,filepath)
% [w,n,comment] = read_nk(filename,scale,filepath)
%
% Read IMD-format nk file. Each text line in the file can be one of the
% following:
%   (1) an empty string
%   (2) a comment line beginning with ';'
%   (3) a line of the form
%         w  n  k
%       wherein
%         w = wavelength (Angstroms)
%         n = real part of refractive index
%         k = imaginary part of refractive index (k >= 0)
%
% Notes:
%
%   A directory of nk files for a variety of materials can be freely
%   downloaded from http://cletus.phys.columbia.edu/~windt/idl/ . (Extract
%   nk.dir from the IMD installation download.)
%
%   The MATLAB function interp1 (linear interpolation) can be used to
%   change the refractive index's wavelength sampling. (interp1 will
%   internally sort the data by wavelength.)
%
% syntax:
%   [w,n,comment] = read_nk(filename,scale,filepath)
%
% input arg's:
%
%   filename: file name; string
%
%   scale: wavelength scaling factor; real (optional; default = 1)
%
%   filepath: file path prefix; string (optional; default = '.\'; filepath
%   must be '\' terminated.)
%
% output arg's
%
%   w: wavelengths, in units of A/scale (e.g., set scale=0.1 for nm, 0.0001
%   for micron); size-[num_w,1] real
%
%   n: complex refractive index (imag(n) >= 0); size-[num_w,1] complex
%
%   comment: file comment lines, size-[1,?] cell array of string
%
% Version 05/12/2005
% Author: Kenneth C. Johnson
% software.kjinnovation.com
% This file is Public Domain software.

if nargin<3
    filepath='.\';
    if nargin<2
        scale=1;
    end
end
fid=fopen([filepath filename]);
if fid==-1
    error(['Could not open file ' filepath filename]);
end
comment={};
wnk=zeros(3,0);
next_line=fgetl(fid);
while ~isequal(next_line,-1)
    if length(next_line)>0
        if isequal(next_line(1),';')
            comment{end+1}=next_line;
        else
            wnk(:,end+1)=sscanf(next_line,'%g');
        end
    end
    next_line=fgetl(fid);
end
w=scale*wnk(1,:).';
n=(wnk(2,:)+i*wnk(3,:)).';
f=find(diff(w)==0);
if ~isempty(f)
    w(f)=[];
    n(f+1)=(n(f)+n(f+1))/2;
    n(f)=[];
end
