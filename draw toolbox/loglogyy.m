function [ hax, h1, h2 ] = loglogyy( x1, y1, x2, y2, varargin )
%[hax, h1, h2] = loglogyy(x1, y1, x2, y2)
%   just like plotyy() but with loglog axis

ha1 = axes;
h1 = loglog(ha1, x1, y1); grid
ha2 = axes;
h2 = loglog(ha2, x2, y2, '-g'); grid
set(ha2,'color','none')
set(ha2,'yaxisLocation','right')
set(ha2,'position',get(ha1,'position'))

hax = [ha1 ha2];


end

