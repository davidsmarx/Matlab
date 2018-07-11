function d = distancepoint2plane(p, plane)
% d = distancepoint2plane(p, plane)
%
% p = point [x y z], where x, y, z can be column vectors
% plane = [a; b; c] where z = a*x + b*y + c = [x y 1] * plane;
%
% see Boas p. 219

[npoints, n3] = size(p);
if n3 ~= 3, error('points must be [x y z]'); end

% normal to plane
n = [-plane(1); -plane(2); 1]./sqrt(plane(1).^2 + plane(2).^2 + 1);

% vector from p to any point on the plane, say 0, 0, c
% pq is 3 x npoints
pq = p - repmat([0 0 plane(3)],npoints,1);

% distance = n . pq, result is a row vector
%d = abs(n' * pq);
d = pq * n;
