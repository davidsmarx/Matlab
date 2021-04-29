function [ hax, h1, h2 ] = plotyygeneral( x1, y1, bLog1, x2, y2, bLog2, varargin )
% [ hax, h1, h2 ] = plotyygeneral( x1, y1, bLog1, x2, y2, bLog2)
%   just like plotyy() but specify each y axis is Log or LInear
% 
% NOTE: need to rewrite this using yyaxis left; yyaxis right

ha1 = axes;
h1 = PlotY(ha1, x1, y1, bLog1); grid on

ha2 = axes;
h2 = PlotY(ha2, x2, y2, bLog2); grid on
set(h2,'color','g')

set(ha2,'color','none')
set(ha2,'yaxisLocation','right')

hax = [ha1 ha2];

end

function h = PlotY(ha, x, y, bLog)

    if bLog,
        h = semilogy(ha, x, y);
    else
        h = plot(ha, x, y);
    end

end