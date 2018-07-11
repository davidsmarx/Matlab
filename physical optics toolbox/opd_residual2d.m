function [linResidual, quadResidual, rms, x, y, linxx, quadxx] =...
    opd_residual2d(varargin)
% [linResidual, quadResidual, rms, x, y, linxx, quadxx] = opd_residual2d(delaylinelist,opd)
% [linResidual, quadResidual, rms, x, y, linxx, quadxx] = opd_residual2d(xfov,yfov,opd)
%
% x, y (rad) are angle in the field of view. y-fov corresponds to delayline
% movement. y-fov -> tilt about x-axis, x-fov -> tilt about y-axis.
%
% linxx = [DC XlinearSlope YlinearSlope]
% quadxx = [DC XlinearSlope YlinearSlope radialQuadCoef]
%
% this routine is obsolete as of 1/13/05. use opd_residual() instead.

% constant parameters:
N = 128; % size of 2-d grid for interpolating opd

% parse inputs and create the 2-d grid
switch nargin,
    case 2,
        % only 1-d delayline movement is input
        % assume opd is constant in x (orthogonal to baseline)
        % and interpolate a 2-d response
        delaylinelist = varargin{1};
        opd = varargin{2};
        
        % convert delaylinelist to field of view angles
        baseline = 10; % (meters)
        fov = asin(delaylinelist./baseline); % (rad)
        % make an evenly spaced grid over the circular FOV
        x = linspace(fov(1),fov(end),N);
        y = x;
        [X, Y] = meshgrid(x,y);
        R = sqrt(X.^2 + Y.^2);

        % use only the fov within a circle
        Ruse = R<=max(fov);

        % interpolate the opd values over the circle
        % movements over the field of view in the y-direction (tilt about
        % x-axis) corresponds to movement in the delayline. So, the OPD map
        % is constant in xfov and = opd(yfov).
        OPD = interp1(fov,opd,Y,'pchip');
        
    case 3,
        % opd is 2-d matrix and x and y articulations are given
        xfov = varargin{1};
        yfov = varargin{2};
        Opd  = varargin{3};

        [X, Y] = meshgrid(xfov,yfov);
        R = sqrt(X.^2 + Y.^2);
        % use only the fov within a circle
        Ruse = R<=max(X(:));
        OPD = Opd;
        x = xfov; y = yfov;
%         % interpolate onto fine grid
%         x = linspace(xfov(1),xfov(end),N);
%         y = linspace(yfov(1),yfov(end),N);
%         [X, Y] = meshgrid(x,y);
%         R = sqrt(X.^2 + Y.^2);
%  
%         % use only the fov within a circle
%         Ruse = R<=max(X(:));
% 
%         % interpolate the opd values over the circle
%         OPD = interp2(xfov,yfov,Opd,X,Y,'cubic',0);
      
    otherwise,
        error('usage mismatch');
end

% fit and rms the residual
linResidual = zeros(size(X));
quadResidual = zeros(size(X));
% [linxx,  linResidual(Ruse),  fval] = lfit_fdepcoeflin( [X(R<Ruse) Y(R<Ruse)],OPD(R<Ruse));
% [quadxx, quadResidual(Ruse), fval] = lfit_fdepcoefQuad([X(R<Ruse) Y(R<Ruse)],OPD(R<Ruse));
[linxx,  linResidual(Ruse),  fval] = lfit_fdepcoeflin( [X(Ruse) Y(Ruse)],OPD(Ruse));
[quadxx, quadResidual(Ruse), fval] = lfit_fdepcoefQuad([X(Ruse) Y(Ruse)],OPD(Ruse));
rms = sqrt(mean(quadResidual(Ruse).^2));

Z = zernikefit(X(Ruse),Y(Ruse),OPD(Ruse),15);
% r = linspace(0,max(fov),N);
% t = [0:99]'*2*pi/N;
% 
% [R, T] = meshgrid(r,t);
% X = R .* cos(T);
% Y = R .* sin(T);
% 
% % interpolate the opd over the whole grid
% OPD = interp1(fov,opd,X,'pchip');
% 
% % calculate linear and quadratic fits
% [linxx, linResidual, fval] = lfit_fdepcoeflin( [X(:) Y(:)],OPD(:));
% [quadxx,quadResidual,fval] = lfit_fdepcoefQuad([X(:) Y(:)],OPD(:));
% 
% rms = sqrt(mean(quadResidual(:).^2))
% rms = sqrt(mean(quadResidual(:).^2.*R(:).*diff(r([1 2])))./mean(R(:).*diff(r([1 2]))))

