function d = distancepoint2line(p, line)
% d = distancepoint2line(p, line)
%
% p = [x; y]
% line = [x1 x2; y1 y2]
%
% see http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/

p1 = line(:,1); p2 = line(:,2);
dp1top2_2 = (p1 - p2)' * (p1 - p2);

if dp1top2_2 <= 0.0,
    % p1 == p2
    d = sqrt( (p - p1)' * (p - p1) );
    
else,
    u = (p(1) - p1(1))*(p2(1) - p1(1)) + (p(2) - p1(2))*(p2(2) - p1(2));
    u = u./dp1top2_2;
    
    % the point where the line and the tangent intersect
    pi = p1 + u*(p2 - p1);
    
    d = sqrt( (pi - p)' * (pi - p) );
    
end