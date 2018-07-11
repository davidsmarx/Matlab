function [U, R, Ro] = z_getglobalcoords(zchan,nsurf,V)
% [U, R, Ro] = z_getglobalcoords(zchan,nsurf,V)
% nsurf = surface to get coordinate system
% V = 3 x N matrix, each column is a direction vector in the local
%     coordinate system of surface nsurf
% U = 3 x N matrix, each column is the direction vector in global
%     coordinate system for the corresponding column in V
% R = the rotation matrix to go from local coords to global coords
% Ro= the offset vector. U = R*V + Ro;

[nr, nc] = size(V);
if nr ~= 3, error('input V must be 3 x N'); end

try
cmdstr = ['GetGlobalMatrix, ' num2str(nsurf)];
RR = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[3 inf]);

R = RR(:,1:3)';
Ro= RR(:,4);

U = R*V + Ro*ones(1,nc);

catch,
   disp('ERROR in z_getglobalcoords!');
   keyboard;
end


return


