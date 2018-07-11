function M = make_cc_mask_movie(node,cc,mask,masktheta,sp)
% M = make_cc_mask_movie(node,cc,mask,masktheta,sp)
%
% node = SIM label of corner cube (1 = without delay line 'ew', 2 = with
% 'ns')
% cc = corner cube struct
% mask = [dx dy x0]; OR
% mask = struct('mask',[xb, dy],'bpoly',[x1 .. xn; y1 .. yn])
% mask = [x1 .. xn; y1 .. yn] polygon vertices
%
% masktheta = mask rotation
% sp = SIM parameter struct

global MM;

switch node,
    case {1,3},
        direction = 'ew';
    case {2,4},
        direction = 'ns';
    otherwise,
        error('unknown node');
end

% nested function for drawing mask
    function hh = DrawMask
        if isequal(class(mask),'struct'),
            hh = drawmask_fsmpoly(gca,mask.mask/MM,mask.bpoly/MM,...
                masktheta,direction);
            
        elseif isequal([3, 1],size(mask)),
            hh = drawmask_rects(gca,mask(1)/MM,mask(2)/MM,mask(3)/MM,...
                masktheta,direction);
        else,
            hh = drawmask_polygons(gca,mask/MM);
        end
    end


% make movie of articulating corner cube
fovt = linspace(0,2*pi,21);
fovr = sp.fov;

hf = figure;

for ii = 1:length(fovt),
    
    fovx = fovr*cos(fovt(ii));
    fovy = fovr*sin(fovt(ii));
    [cc.rotmatrix, cc.edgedir] = corner_cube_orientations(node,fovx,fovy);
    
    figure(hf);
    hc = draw_CCube_projection(gca,cc,'FaceAlpha',0.25);
    hold on
    hh = DrawMask;
    axis image
    hold off
    grid on
    
    set(gca,'xlim',[-10 10])
    set(gca,'ylim',[-10 10])
    set(gca,'zlim',[0 15])
    
    pause;
    M(ii) = getframe;
end

movie(M)

[filename, pathname, fi] = uiputfile({'*.avi'},'Enter Filename for Movie');
if ~isequal(filename,0),
    movie2avi(M,[pathname filename],'fps',3,'compression','none');
end

end