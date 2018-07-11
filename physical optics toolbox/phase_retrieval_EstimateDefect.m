function [Arefguess, Zrefguess, maxdepth, radius] = phase_retrieval_EstimateDefect(...
    Adiff, dx, dy, zprop, wavelength, sOptions)
% [Arefguessinit, maxdepth, radius] = CreateArefInit(...
%    Adiff, dx, dy, zprop, wavelength, sOptions)
%
% Adiff is the complex field, but only amplitude is used
%
% sOptions struct:
%     'debug' = 'off' (default), 'on'
%
% Arefguess = complex field with amplitude = 1. You should multiply
%      this result by the measured reference field amplitude
% Zrefguess = estimated object height profile.
%      Arefguess = exp(-j*k*2*Zrefguess)

persistent hfigInitDiameter;

if exist('sOptions','var'), sOptions = CheckOptions(sOptions); else sOptions = CheckOptions; end
    function sOptout = CheckOptions(sOptin)
        sOptout = struct('debug','off');
        if nargin >= 1,
            sInfields = fieldnames(sOptin);
            for ifield = 1:length(sInfields),
                sOptout.(sInfields{ifield}) = sOptin.(sInfields{ifield});
            end
        end
    end % CheckOptions

% the physical grid
[Ny, Nx] = size(Adiff);
x = (-Nx/2:Nx/2-1)*dx;
y = (-Ny/2:Ny/2-1)'*dy;
[X, Y] = meshgrid(x,y);
R = sqrt(X.^2 + Y.^2);

% estimate the curvature and diameter of the bowl to create initial guess
roc = 2*zprop;

% find dark pixels and bright pixels.
% Create bw image from those pixels.
% Use morphologic closing and choose the blob near the middle.
% Trace the edge of the blob use volume to estimate area of defect
meanAdiff = mean(Adiff(:));

% use only the central half of the image because we are given that the
% defect is at the center of the image
icenter = R < max(x)/4;

%         epsilon = .01*meanAdiff;
%         idark = Adiff < meanAdiff - epsilon & abs(X) < 0.25*Nx*dx & abs(Y) < 0.25*Ny*dy;
% dark pixels:
minAdiff = min(Adiff(icenter));
idark = Adiff < [0.2 0.8]*[meanAdiff minAdiff]' & R < max(x)/2;

% bright pixels:
maxAdiff = max(Adiff(icenter));
ibright = Adiff > [0.5 0.5]*[meanAdiff maxAdiff]' & R < max(x)/2;

% plot the
Abw = zeros(Ny, Nx, 3); % true rgb image
Abw(:,:,3) = idark; % blue
Abw(:,:,1) = ibright; % red

% combine bright and dark
% parameters are for 128x128 image
%Bbw = imclose(idark | ibright, strel('disk',10));
Bbw = imclose(idark, strel('disk',10));
Bbw = imopen(Bbw, strel('disk',1*Nx/128));

% find a circle that encloses the center blob
[Blabel, numobjects] = bwlabel(Bbw,8);
% if there is only one object, use it
switch numobjects
    case 0
        % do the old method
        % old method:
        darkarea_pixels = sum(idark(:)); % # of dark pixels
        darkradius = sqrt(dx*dy*darkarea_pixels / pi);
        
        % no centroid information
        xcentroid = sum(X(idark))./darkarea_pixels;
        ycentroid = sum(Y(idark))./darkarea_pixels;
        
    case 1
        % find circle that encloses the object
        iobjuse = 1;
        numpix = sum(Blabel(:) == iobjuse);
        xcentroid = sum(X(Blabel(:) == iobjuse))./numpix;
        ycentroid = sum(Y(Blabel(:) == iobjuse))./numpix;

        % calculate radius of the found object
        Robj = sqrt( (X(Blabel == iobjuse) - xcentroid(iobjuse)).^2 + ...
            (Y(Blabel == iobjuse) - ycentroid(iobjuse)).^2 );
        darkradius = max(Robj);

        % calculate centroid of the found object
        % find the center of the defect --- weighted centroid inside the
        % darkradius
        % [xdef, ydef, Adef] = filterdata(R < darkradius, X, Y, Adiff);
        % xcentroid = sum(xdef.*Adef)./sum(Adef);
        % ycentroid = sum(ydef.*Adef)./sum(Adef);
        xcentroid = xcentroid(iobjuse);
        ycentroid = ycentroid(iobjuse);
    
    otherwise
        % find the largest object whose centroid is closest to image center
        [numpix, xcentroid, ycentroid] = deal(zeros(numobjects,1));
        for iobj = 1:numobjects
            numpix(iobj) = sum(Blabel(:) == iobj);
            xcentroid(iobj) = sum(X(Blabel(:) == iobj))./numpix(iobj);
            ycentroid(iobj) = sum(Y(Blabel(:) == iobj))./numpix(iobj);
            %rcentroid_2(iobj) = [xcentroid ycentroid]*[xcentroid ycentroid]';
        end
        rcentroid_2 = [xcentroid ycentroid] * [xcentroid ycentroid]';
        maxsize = max(numpix);
        iobjmax = find(numpix == maxsize);
        if length(iobjmax) > 1,
            [minr, iobjmin] = min(rcentroid_2(iobjmax));
            iobjuse = iobjmax(iobjmin);
        else
            iobjuse = iobjmax;
        end
        
        % calculate radius of the found object
        Robj = sqrt( (X(Blabel == iobjuse) - xcentroid(iobjuse)).^2 + ...
            (Y(Blabel == iobjuse) - ycentroid(iobjuse)).^2 );
        darkradius = max(Robj);
        
        % calculate centroid of the found object
        % find the center of the defect --- weighted centroid inside the
        % darkradius
        % [xdef, ydef, Adef] = filterdata(R < darkradius, X, Y, Adiff);
        % xcentroid = sum(xdef.*Adef)./sum(Adef);
        % ycentroid = sum(ydef.*Adef)./sum(Adef);
        xcentroid = xcentroid(iobjuse);
        ycentroid = ycentroid(iobjuse);

end % switch numobjects

if strcmp(sOptions.debug, 'on'), Debug; end
    function Debug        
        if isempty(hfigInitDiameter)
            hfigInitDiameter = tamarfigure;
        end
        figure(hfigInitDiameter),
        
        subplot(3,1,1)
        imagesc(x, y, Adiff), colormap gray, axis image
        hold on
        ttmp = linspace(-pi,pi);
        plot(xcentroid + darkradius*cos(ttmp), ycentroid + darkradius*sin(ttmp), 'r')
        hold off
        
        subplot(3,1,2)
        imshow(Abw) %, subplot(2,1,2), imshow(Bbw)

        subplot(3,1,3)
        imshow(Bbw)
    end % function Debug

% so defect is roughly darkradius with roc
Zrefguess = zeros(Ny,Nx);
Roff = sqrt( (X-xcentroid).^2 + (Y-ycentroid).^2);

% depth < 0 ==> divot, which is the assumption here
% raised cosine scaled so that Taylor series quadratic term has the
% appropriate radius of curvature
% old simple quadratic:
%depth = -darkradius.^2./(2*roc);
%Zinit(Roff < darkradius) = -depth + Roff(Roff < darkradius).^2./(2*roc);
if 0
    % raised cosine:
    depth = -2*darkradius.^2./(pi^2*roc);
    idef = Roff < darkradius;
    Zrefguess(idef) = depth*0.5*(1 + cos(pi*Roff(idef)./darkradius));
else
    % gaussian with darkradius = 1/e^2 point
    depth = -0.25*darkradius.^2./roc;
    alpha_2 = 2*depth*roc;
    Zrefguess = depth*exp(-Roff.^2./abs(alpha_2));
end

% figure, mysurf(X, Y, Zrefguess);
% title(['depth = ' num2str(depth,'%e') ', \alpha^2 = ' num2str(alpha_2,'%e')]);

Arefguess = exp(j*2*(2*pi./wavelength)*Zrefguess);

maxdepth = depth;
radius = darkradius;

end % CreateArefInit
