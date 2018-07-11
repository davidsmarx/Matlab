function [linResidual, quadResidual, quadRms, x, y, linxx, quadxx, Z] =...
    opd_residual(opd,simparms,varargin)
% [linResidual, quadResidual, quadRms, x, y, linxx, quadxx, Z] = ...
%     opd_residual(opd,simparms,options)
%
% opd is a matrix of opd values defined on the field of view grid defined
% in simparms.
% or opd is a vector corresponding to the y-fov, in which case the opd is
% replicated along the x-fov (assumes corner cube tilt has no effect) to
% create the 2-d grid.
%
% sim_parms is a struct with the following fields and defaults:
%    baseline = SIM baseline (default = 10m)
%    fov      = the angular radius of the FOR (default = 1degree--narrow
%               angle)
%    fovX0    = center of FOR in x-direction (default = 0)
%    fovY0    = center of FOR in y-direction (default = 0)
%    nfov     = number of grid points in one dimension (default = 21, total
%               number of grid points = nfov x nfov, except that only those
%               within the fov radius are used.
% sim_parms = calc_2dopd_ce_hetpow; will return a struct with the default
% values.
%
% x, y (rad) are angle in the field of view. y-fov corresponds to delayline
% movement. y-fov -> tilt about x-axis, x-fov -> tilt about y-axis.
%
% linxx = [DC XlinearSlope YlinearSlope]
% quadxx = [DC XlinearSlope YlinearSlope radialQuadCoef]
%
% Z = first 28 Zernike coefficients of OPD map
%
% options:


% % only 1-d delayline movement is input
%         % assume opd is constant in x (orthogonal to baseline)
%         % and interpolate a 2-d response
%         delaylinelist = varargin{1};
%         opd = varargin{2};
%         
%         % convert delaylinelist to field of view angles
%         baseline = 10; % (meters)
%         fov = asin(delaylinelist./baseline); % (rad)
%         % make an evenly spaced grid over the circular FOV
%         x = linspace(fov(1),fov(end),N);
%         y = x;
%         [X, Y] = meshgrid(x,y);
%         R = sqrt(X.^2 + Y.^2);
% 
%         % use only the fov within a circle
%         Ruse = R<=max(fov);
% 
%         % interpolate the opd values over the circle
%         % movements over the field of view in the y-direction (tilt about
%         % x-axis) corresponds to movement in the delayline. So, the OPD map
%         % is constant in xfov and = opd(yfov).
%         OPD = interp1(fov,opd,Y,'pchip');
        
% opd is 2-d matrix on grid defined by simparms
tmpfov = linspace(-simparms.fov,simparms.fov,simparms.nfov)';
xfov = tmpfov + simparms.fovX0;
yfov = tmpfov + simparms.fovY0;

% if opd is 1-d, replicate it in x-direction
if isvector(opd),
    opd = opd(:) * ones(1,length(xfov));
end

[X, Y] = meshgrid(xfov,yfov);
R = sqrt((X-simparms.fovX0).^2 + (Y-simparms.fovY0).^2);
% use only the fov within a circle
Ruse = R<=simparms.fov;
x = xfov; y = yfov;

% fit and rms the residual
linResidual = zeros(size(X));
quadResidual = zeros(size(X));
% [linxx,  linResidual(Ruse),  fval] = lfit_fdepcoeflin( [X(R<Ruse) Y(R<Ruse)],OPD(R<Ruse));
% [quadxx, quadResidual(Ruse), fval] = lfit_fdepcoefQuad([X(R<Ruse) Y(R<Ruse)],OPD(R<Ruse));
[linxx,  linResidual(Ruse),  fval] = lfit_fdepcoeflin( [X(Ruse) Y(Ruse)],opd(Ruse));
[quadxx, quadResidual(Ruse), fval] = lfit_fdepcoefQuad([X(Ruse) Y(Ruse)],opd(Ruse));
quadRms = sqrt(mean(quadResidual(Ruse).^2));

opduse = ~isnan(opd);
rmsresidual = zeros(1,28);
for ii = 1:28,
    [Z, restmp, rmstmp, R] = zernikefit(X,Y,opd,ii,simparms.fov);
    rmsresidual(ii) = rmstmp;
end

return
