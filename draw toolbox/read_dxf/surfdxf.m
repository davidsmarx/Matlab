% FUNCTION hh = surfdxf(glt)
%
% This function draws a surface plot
% of a file read in with 'read_dxf'
%
function hh = surfdxf(glt,varargin)
  
  s = size(glt);
  gl2 = reshape(glt,s(1),s(2)*s(3));
  x = gl2(1,:);
  y = gl2(2,:);
  z = gl2(3,:);
  
  tri = 1:(s(2)*s(3));
  tri = reshape(tri,s(2),s(3))';
  
  if(nargin > 1)
    h = trisurf(tri,x,y,z,varargin);
  else
    h = trisurf(tri,x,y,z);
  end
  
  if(nargout==1) 
    hh=h;
  end
  
