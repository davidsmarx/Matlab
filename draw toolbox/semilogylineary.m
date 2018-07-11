function [ hax, h1, h2 ] = semilogyy( x1, y1, x2, y2, varargin )
%[hax, h1, h2] = semilogyy(x1, y1, x2, y2)
%   just like plotyy() but with semilogy axis

ha1 = axes;
h1 = semilogy(ha1, x1, y1); grid
ha2 = axes;
h2 = semilogy(ha2, x2, y2, '-g'); grid
set(ha2,'color','none')
set(ha2,'yaxisLocation','right')

hax = [ha1 ha2];


end

