function bIn = InTriangle(pt, p1, p2, p3)
% bIn = InTriangle(pt, p1, p2, p3)
%
% test if pt is inside the triangle formed by p1, p2, p3
%
% p = [x; y] or p = [x; y; z] or p = [x; y; 0]

if length(pt) == 2, pt = [pt; 0]; end
if length(p1) == 2, p1 = [p1; 0]; end
if length(p2) == 2, p2 = [p2; 0]; end
if length(p3) == 2, p3 = [p3; 0]; end

% shift all points to place pt at the origin
p1t = p1 - pt;
p2t = p2 - pt;
p3t = p3 - pt;

% use cross products, if all the same sign, then point is inside triangle
cp1 = cross(p1t, p2t);
cp2 = cross(p2t, p3t);
cp3 = cross(p3t, p1t);

bIn = false;
if all([cp1(3) cp2(3) cp3(3)] > 0), bIn = true; end
if all([cp1(3) cp2(3) cp3(3)] < 0), bIn = true; end

