function phiref = phase_retrieval_HIO(...
    Iref, Idiff, dx, dy, zprop, wavelength, sOptions) % Arefguessinit)
% [zsurf, phiref] = phase_retrieval_HIO(
%        Iref, Idiff, x, y, zprop, wavlength, options)
%
% Iref, Idiff are intensity images ( = amp^2 )
% dx, dy are grid spacing
% zprop = propagation distance from 'ref' plane to 'diff' plane
% sOptions = struct with fields (all optional)
%     'ArefInit' = initial guess at complex field in the ref plane (default
%         is empty)
%     'AdiffInit' = initial guess at complex field in the diff plane
%         (default is empty)
%     'debug' (default = 'off')
%
% reference:
%    J.R. Fienup, "Phase Retrieval Algorithms: A Comparison," Applied
%    Optics, vol. 21, p. 2758, 1982.

%unitsdefinitions;
%global MM UM;

% constants, maybe should be options
LOW_FREQ_CUTOFF = 0.25;
BETA = 0.5; %

% validate inputs
[Adiff, Aref, Nx, Ny] = ValidateInputs(Iref, Idiff, dx, dy);

if ~exist('sOptions','var'), sOptions = struct; end
sOptions = CheckOptions(sOptions);

% if the total energy in the ref field and diffracted field are very
% different, renormalize
%if abs(sum(Iref(:)) - sum(Idiff(:))) > 0.1*sum(Iref(:)),
    Adiff = sqrt(sum(Iref(:))./sum(Idiff(:)))*Adiff;
    %warning('renormalizing Adiff');
%end

% the physical grid
x = (-Nx/2:Nx/2-1)*dx;
y = (-Ny/2:Ny/2-1)'*dy;
% [X, Y] = meshgrid(x,y);
% R = sqrt(X.^2 + Y.^2);
% T = atan2(Y,X);

% check for an initial guess, otherwise need to create an initial field
% need an initial guess of the diffracted field, so if the input initial
% estimate is the reference field, propagate it to the diffracted field
if isempty(sOptions.ArefInit) && isempty(sOptions.AdiffInit),
    % estimate the object phase, then multiply by measured object amplitude
    Arefguessinit = phase_retrieval_EstimateDefect(...
        Adiff, dx, dy, zprop, wavelength, sOptions);
    Arefguessinit = Aref.*Arefguessinit;
    [Adiffguess, curvout, dxout, dyout] = propagate(Arefguessinit, dx, dy, 0, wavelength, zprop, dx, dy);
elseif ~isempty(sOptions.ArefInit)
    [Adiffguess, curvout, dxout, dyout] = propagate(sOptions.ArefInit, dx, dy, 0, wavelength, zprop, dx, dy);
else % must be sOptions.AdiffInit is not empty
    Adiffguess = sOptions.AdiffInit;
end

% iterate
persistent hfigIterate;
if strcmp(sOptions.debug, 'on'),
    if isempty(hfigIterate)
        hfigIterate = figure;
    end
end
notdone = true;
fom = zeros(100,1);
icnt = 0;
while notdone

    % propagate back to object plane
    %Adiffguess = SpatialFilter(Adiffguess);
    [Arefprop, curvout, dxouttmp, dyouttmp] = propagate(Adiffguess, dx, dy, 0, wavelength, -zprop, dx, dy);

    % make correction to object field
    % gp_k = Arefguess
    % g_k = Arefguess_last = Aref (amplitude)
    % g_k+1 = next Arefguess
    % g_k+1 = gp_k where constraints are satisfied
    %         g_k - beta*gp_k where contraints are violated
    amp_p_k = abs(Arefprop);
    amp_kp1 = Aref - BETA*(amp_p_k - Aref);
    amp_kp1(amp_kp1 < 0) = 0;
    Arefguess = amp_kp1 .* (Arefprop./amp_p_k);
    
    Arefguess = ApplySpatialFreqCutoff(Arefguess, LOW_FREQ_CUTOFF);
    
    % propagate to diffracted plane
    [Adiffguess, curvout, dxout, dyout] = propagate(Arefguess, dx, dy, 0, wavelength, zprop, dx, dy);

    % measured image is low pass filtered image of diffracted field
    %Adiffguessimage = ApplySpatialFreqCutoff(Adiffguess, 0.0625);
    
    if strcmp(sOptions.debug,'on'),
        figure(hfigIterate),
        set(gcf,'position',[560   251   560   697]),
        subplot(3,1,1), imagesc(x, y, abs(Arefprop)  - Aref), colorbar
        title(['rms error = ' num2str(rms(abs(Arefprop(:))-abs(Aref(:))),'%.3f')])
        ppp = angle(Arefprop(:,Nx/2+1));
        zpp = unwrap(ppp)/(0.6*4*pi/wavelength);
        subplot(3,1,2), plot(x,zpp), grid
        title(['iteration #' num2str(icnt)])
        subplot(3,1,3), 
        %imagesc(x, y, abs(Adiffguess) - Adiff), colorbar
        %plot(y, [abs(Adiffguess(:,Nx/2+1)) abs(Adiff(:,Nx/2+1)) abs(Adiffguessimage(:,Nx/2+1))]), grid, legend('guess','input','image')
        plot(y, [abs(Adiffguess(:,Nx/2+1)) abs(Adiff(:,Nx/2+1))]), grid, legend('guess','input')
        title(['rms error = ' num2str(rms(abs(Adiffguess(:))-abs(Adiff(:))),'%.3f')])
        pause(0.5);
    end
    
    % make correction to amplitude of diffracted field
    % gp_k = Adiffguess
    % g_k = Adiffguess_last
    % g_k+1 = gp_k where constraints are satisfied
    %         g_k - beta*gp_k where contraints are violated
    amp_p_k = abs(Adiffguess);
    %amp_p_k = abs(Adiffguessimage);
    amp_kp1 = Adiff - BETA*(amp_p_k - Adiff);
    amp_kp1(amp_kp1 < 0) = 0;
    Adiffguess = amp_kp1 .* (Adiffguess./amp_p_k);

    Adiffguess = ApplySpatialFreqCutoff(Adiffguess, LOW_FREQ_CUTOFF);

    % increment the iteration counter
    icnt = icnt + 1;

    % exit criteria (needs work)
    fom(icnt) = max(abs( abs(Adiffguess(:)) - Adiff(:) ));
    if icnt >= 10, notdone = false; end
    
end

if strcmp(sOptions.debug,'on')
    save C:\PhaseRetrievalHIOlog.mat;
end

% return values:
phiref = angle(Arefguess);


end % main

function [Adiff, Aref, Nx, Ny] = ValidateInputs(Iref, Idiff, dx, dy)

if any(size(Iref) ~= size(Idiff)),
    error('input arrays must be same size');
end
[Ny, Nx] = size(Iref);

if dx ~= dy,
    error('input grid spacing must be equal');
end

Adiff = sqrt(Idiff);
Aref = sqrt(Iref);

end % ValidateInputs

function sOptions = CheckOptions(sIn)
% sOptions = struct with fields (all optional)
%     'ArefInit' = initial guess at complex field in the ref plane (default
%         is empty)
%     'AdiffInit' = initial guess at complex field in the diff plane
%         (default is empty)
%     'debug' (default = 'off')

% default values:
sOptions = struct(...
    'ArefInit',[],...
    'AdiffInit',[],...
    'LibPropParms',struct(),...
    'debug','off');

% copy any input options to the options struct
sInfields = fieldnames(sIn);
for ifield = 1:length(sInfields),
    sOptions.(sInfields{ifield}) = sIn.(sInfields{ifield});
end

end % CheckOptions

function Aout = ApplySpatialFreqCutoff(Ain, Fc_normalized)
    % since the images are resampled, expanded versions of the input, the
    % answer must be within spatial frequency bounds
    [Ny, Nx] = size(Ain);
    Afft = fftshift(fft2(Ain));
    
    Fclow = @(N) floor(N/2 - N/2*Fc_normalized);
    Fchigh = @(N) ceil(N/2 + N/2*Fc_normalized);

    %     Fclow = Ny/2 - Ny/2*Fc_normalized;
    %     Fchigh = Ny/2 + Ny/2*Fc_normalized;
    %     Afft([1:Ny/4 3*Ny/4+1:end],:) = 0;

    % raised cosine to soften edges and reduce ringing
    Ftaperwidth = 8; % pixels
    % separate N pixels into inside and outside the passband
    iblockneg = 1:Fclow(Nx)-1;
    ipass = Fclow(Nx):Fchigh(Nx);
    iblockpos = Fchigh(Nx)+1:Nx;
    Npass = length(ipass);
    alphax = Ftaperwidth/(Npass/2);
    alphay = Ftaperwidth/(Npass/2);
    winx = [
        zeros(length(iblockneg),1);
        tukeywin(length(ipass),alphax);
        zeros(length(iblockpos),1)
        ];
    winy = winx; %[zeros(Fclow(Ny)-Ftaperwidth,1); tukeywin(Ny-Fclow(Ny)-(Ny-Fchigh(Ny))+2*Ftaperwidth,alphay); zeros((Ny-Fchigh(Ny))-Ftaperwidth,1)];
    [WX, WY] = meshgrid(winx, winy);
    W = WX.*WY;
    
    %     Afft([1:Fclow(Ny) Fchigh(Ny)+1:end],:) = 0;
    %     Afft(:,[1:Fclow(Nx)-1 Fchigh(Nx):end]) = 0;
    Afft = W .* Afft;

    Aout = ifft2(fftshift(Afft));
    
end
