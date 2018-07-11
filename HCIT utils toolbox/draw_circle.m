
function hout = draw_circle(ctr, D, lw, col)
% draw_circle(ctr, D, lw, col)

R       = D/2;
h=rectangle('Position', ...
    [ctr-R*[1 1] [D D]],...
    'EdgeColor', col, 'LineWidth',lw,'Curvature',[1 1]);
                    
if nargout > 0,
    hout = h;
end

return
