function [xx, yy, ww] = GaussQuad_CircleDomain(R, Nrings, Nspokes)
% [xx, yy, ww] = GaussQuad_CircleDomain(R, Nrings, Nspokes)

    [x, wf] = gausswts(0, 1, Nrings);
    theta = (0:Nspokes-1)*2*pi/Nspokes;
    
    rr = R.*sqrt(x);
    wf = (2*pi/Nspokes) .* (0.5.*R.^2.*wf);
    
    xx = rr(:) * cos(theta(:)');
    yy = rr(:) * sin(theta(:)');
    ww = wf(:) * ones(1,Nspokes);
    
    xx = xx(:)';
    yy = yy(:)';
    ww = ww(:)';
    
end % GaussQuad

% % example calculating the area of a circle:
% radius = 3;
% [xx, yy, ww] = GaussQuad_CircleDomain(radius, 3, 2)
% area = sum(ww.*ones(size(xx)))
% error =  pi*radius^2 - area