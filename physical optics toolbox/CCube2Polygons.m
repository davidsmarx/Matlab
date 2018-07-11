function [vp1 vp2 vp3 vpboundary] = CCube2Polygons(cc)
% [vp1 vp2 vp3 vpboundary] = CCube2Polygons(cc)
%
% vp1, vp2, vp3 are polygons (nx3 matrix) for each face
% vpboundary is a polygon (nx3 matrix) for the outer boundary
%
% The fields in the Ccube structure are:
%    size
%    shape ('square' (default), 'triangle', 'sim')
%    edgelen(1:3) = relative length of each roofline. (default = [1 1 1]).
%         The length of each roofline will be size*edgelen(ii).
%    xc  = vertex offset in x (default = 0)
%    yc  = vertex offset in y (default = 0)
%    spin = clocking of corner cube (default = 0)
%    gapwidth(1:3) = size of gap along each roofline (default = 0)
%    dihedral(1:3) = dihedral (rad) error of each roofline
%    rotmatrix(1:3,1:3) = specifies tip/tilt (default = identity matrix)
%    edgedir(1:3,1:3) = matrix of three column vectors, where each vector is a
%       direction vector indicating the direction of a roof line.

% make polygons for each face
switch cc.shape,
    case 'sim',
        % first face is the baseplate
        edgetheta = acos(cc.edgedir(3,1));
        edgephi   = atan2(cc.edgedir(2,1),cc.edgedir(1,1));
        NP = 10;
        psi = linspace(0,pi/2,NP+1);
        ctmp = zeros(3,NP);
        for ii = 1:NP,
            ctmp(1:3,ii) =...
                make_rotmatrix(0,0,edgephi)*make_rotmatrix(0,edgetheta,0)*...
                make_rotmatrix(0,0,psi(ii))*...
                make_rotmatrix(0,-edgetheta,0)*make_rotmatrix(0,0,-edgephi)*...
                cc.edgedir(:,2);
        end
        vp1 = [
            0 0 0;
            cc.edgedir(:,3)'*cc.size*cc.edgelen(3);
            flipud(ctmp'*cc.size*cc.edgelen(2));
            %cc.edgedir(:,2)'*cc.size*cc.edgelen(2); % repeats last point
            0 0 0;
            ];
        % the other two faces are rectangles
        vp2 = [
            0 0 0;
            cc.edgedir(:,2)'*cc.size*cc.edgelen(2);
            (cc.edgedir(:,2)'*cc.edgelen(2) + cc.edgedir(:,1)'*cc.edgelen(1))*cc.size;
            cc.edgedir(:,1)'*cc.size*cc.edgelen(1);
            0 0 0;
            ];
        vp3 = [
            0 0 0;
            cc.edgedir(:,1)'*cc.size*cc.edgelen(1);
            (cc.edgedir(:,3)'*cc.edgelen(3) + cc.edgedir(:,1)'*cc.edgelen(1))*cc.size;
            cc.edgedir(:,3)'*cc.size*cc.edgelen(3);
            0 0 0;
            ];
            
        vpboundary = [
            vp1(2:end-2,:);
            vp2(2:end-1,:);
            vp3(2:end-1,:);
            ];
        
    otherwise,
        error(['shape ' cc.shape ' is not yet implemented']);
end

% translate the corner cube according to the vertex offset
vp1 = vp1 + ones(size(vp1,1),1)*[cc.xc cc.yc 0];
vp2 = vp2 + ones(size(vp2,1),1)*[cc.xc cc.yc 0];
vp3 = vp3 + ones(size(vp3,1),1)*[cc.xc cc.yc 0];
