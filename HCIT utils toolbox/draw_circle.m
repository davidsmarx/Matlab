
function hout = draw_circle(ctr, D, lw, col, varargin)
% hout = draw_circle(ctr, D, lw, col, options)
% options are property, value pairs as for rectangle()

R       = D/2;
h=rectangle('Position', ...
    [ctr-R*[1 1] [D D]],...
    'EdgeColor', col, 'LineWidth',lw,'Curvature',[1 1], varargin{:});
                    
if nargout > 0,
    hout = h;
end

return
