% FUNCTION dxf = read_dxf(filename)
% 
% Author:  Steven Michael (smichael@ll.mit.edu)
%
% Date:    3/10/2005
%
% Description
%
%   The following compiled MATLAB file reads in an 
%   ascii DXF file and loads the information into
%   the "dxf" variable
%
% Inputs:
%
%   Filename   :    The filename of the ASCII DXF File
%
% Outputs:
%
%   dxf        :    A 3-D variable of size (NX3X3).  The first 
%                   index N is the number of facets.  The second
%                   index references the 3 vertices of each 
%                   facet, and the third index is the 
%                   (x,y,z) location of the vertex.
%
% Note: 
%
%   This may not work with all DXF files.  The DXF file format
%   is complicated and continuously evolving.  This has worked
%   with every file I happened to use it on, but that will
%   not be true for everyone.
%
%   Also, the function does not load any color or texture
%   information.
%
