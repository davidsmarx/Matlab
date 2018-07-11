function [z, mu, sig] = min_spot_centroid(rcos,rpos)
% [z, mu, sig] = min_spot_centroid(rcos,rpos)
%
% rcos = 3 x N matrix, each column is a direction vector,
% [cos(tx);cos(ty);cos(tz)] returned from a raytrace
% rpos = 3 x N matrix, each column is a position vector, [x;y;z], returned
% from a raytrace. A column of rcos and the corresponding column of rpos
% defines a ray so that the 3-d path of a ray is m*rcos + rpos, for scalar
% m.
%
% z = z-value (optical axis) where the rays minimize the rms difference
% between each ray and the centroid of all the rays.
% mu = centroid at z
% sig = rms of the distance of each ray from the centroid at z
%
% see p. 110 of blue Tamar notebook 5/15/03 - ??
% 12/22/04

N = ValidateInputs(rcos,rpos); % N = # of input rays

% scale direction vector so that z component = 1. Then z*direction vector
% puts the ray at the same axial (z) location for each ray.
rcosp = rcos ./ repmat(rcos(3,:),3,1);

% centroid is mean of direction vectors and mean of position vectors
mucos = mean(rcosp')'; % 3 x 1 column vector
mupos = mean(rpos')';  % 3 x 1 column vector

alp = rcosp - mucos*ones(1,N);
bet = rpos  - mupos*ones(1,N);

z = -sum(diag(alp'*bet)) ./ sum(diag(alp'*alp));
% get centroid at z
mu = z*mucos + mupos;
% difference between rays and centroid
v = z*alp + bet;
sig = sqrt(mean(sum(v.^2)));

%%
% zz = linspace(z-1,z+1,11);
% for ii = 1:11,
%     vv = zz(ii)*alp + bet;
%     sig2(ii) = sqrt(mean(sum(vv.^2)));
% end
% figure, plot(zz,sig2,zz,sigg), grid

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N = ValidateInputs(rcos,rpos)

if any(size(rcos) ~= size(rpos)),
    error('inputs rcos and rpos must have same number of rays');
end

[n3, N] = size(rcos);
if n3 ~= 3,
    error('input matrices must have three rows');
end

% for now assume that each rpos has same z, i.e. each ray starts at same
% axial surface.
if max(abs(rpos(3,:) - mean(rpos(3,:)))) > 1e-6,
    error('each ray must start at same z position');
end

