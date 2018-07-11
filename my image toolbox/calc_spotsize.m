function [radius] = calc_spotsize(x,y,a,lev)
% [radius] = calc_spotsize(x,y,a,lev)
% 
% x, y = vectors defining the x and y grid axes
% a = image
%
% lev = linear intensity threshold [0..1] to determine spot size
% radius is a vector, size(radius) = size(lev)
%
% calc_spotsize normalizes the image a so that its peak = 1, then
% uses contourc() to get the contours at levels lev. The return values
% are the mean radius of each contour

if any(size(a) ~= [length(y) length(x)]),
    error('length(x) = #columns(a) and length(y) = #rows(a)');
end

if length(lev) == 1, clev = [lev lev]; else, clev = lev; end

[nr, nc] = size(a);
    
    [fm, xm, ym] = findpeak2(a);
    a = a./fm;
    
    xma = interp1([1:nc],x,xm);
    yma = interp1([1:nr],y,ym);
    
    C = contourc(x-xma,y-yma,a,[clev]);
    %        C = [level1 x1 x2 x3 ... level2 x2 x2 x3 ...;
    %              pairs1 y1 y2 y3 ... pairs2 y2 y2 y3 ...]

    nc = 1;
    for ii = 1:length(lev),
        levii = C(1,nc);
        npii = C(2,nc);
        x1 = C(1,nc+1:nc+npii);
        y1 = C(2,nc+1:nc+npii);
        radius(ii) = mean(sqrt(x1.^2 + y1.^2));

        nc = nc + npii + 1;

    end
    
return    