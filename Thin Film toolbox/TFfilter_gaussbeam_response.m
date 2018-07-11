function [eta0, etam] =...
   TFfilter_gaussbeam_response(beamdiam,lamlist,n,d,options)
% [eta0, etam, lamlist] = ...
%       TFfilter_gaussbeam_response(beamdiam,wavelengths,n,d,TFparms)
%
%   aoi0 = nominal aoi to filter [rad]
%   beamdiam = beam waist diameter [m]
%   wavelength = list of wavelengths [m]
%   n = n(1) = index of substrate, n(end) = index of incident medium
%   d = thickness of each layer, length(d) = length(n)-2
%   options = struct:
%       'pol' (optional) = 'P' or 'S' (default = 'P')
%       'aoi' (optional) = angle of incidence, (default = 0)
%       'Naoi' (optional) = # of angles for planewave decomposition
%            (default = 400)
%       'R_defocus' (optional) = radius of curvature for defocus
%       'Z_propagate' (optional) = additional propagation (defocus?)

global UM NM P;

% [um] beam waist radius
wx = beamdiam/2;  wy = beamdiam/2; 

% 
Nlam = length(lamlist);
lam0 = mean([min(lamlist) max(lamlist)]); % center wavelength

[pol, aoi0, Naoi, R_defocus, Z_propagate] = CheckOptions;

% AOI range for calculating filter response
un = 6*lam0/(sqrt(2)*pi*wx); % max plane wave for calc filter
aoilist = linspace(aoi0-un,aoi0+un,Naoi); % aoi
if ~isreal(aoilist), warning('error aoilist'); keyboard; end

% for each aoi, get the response
[RD, TD] = TFfilter_response(n, d, aoilist, lamlist, pol);

% figure, imagesc(aoilist/P,lamlist/NM,abs(RD)),
% xlabel12('AOI (deg)'), ylabel12('Wavelength (nm)')
% title('Reflectance Amplitude','fontsize',14)

% for each wavelength, use angular spectrum to calculate coupling coeff.

% filter normal
naoi = [sin(aoi0); 0; cos(aoi0)];

% create grid in angular spectrum domain
% Nsig = 24; % umax: S(umax) = exp(-(Nsig^2)/2)
Nu   = 1024;                   Nv   = 1024;
%umax = Nsig/(sqrt(2)*pi*wx);  vmax = Nsig/(sqrt(2)*pi*wy);
%umax = 1./(sqrt(2).*max(lamlist)); vmax = umax;
epsilon = 1e-8;
umax = sqrt(-0.5*log(epsilon))./(pi*wx); vmax = sqrt(-0.5*log(epsilon))./(pi*wy);
du   = 2*umax/Nu;             dv   = 2*vmax/Nv;
u    = [-Nu/2:Nu/2-1]'*du;    v    = [-Nv/2:Nv/2-1]'*dv;
[U, V] = meshgrid(u,v);

% calculate plane wave spectrum of gaussian beam
S = pi*sqrt(wx*wy)*exp(-(pi.^2)*((wx.*U).^2 + (wy.*V).^2));
SS0 = abs(sum(sum(S.^2))).^2;

% figure, imagesc(u, v, abs(S))
% xlabel12('U'), ylabel12('V'), colorbar, title('|S|','fontsize',12)

[eta0, etam] = deal(zeros(Nlam,1));
for nl = 1:Nlam,
   lam  = lamlist(nl);

   % estimate filter response for each plane wave
   % propagation vector, u = sin(theta)/lam:
   kz = sqrt((1/lam).^2-U.^2-V.^2);
   kk = lam*[U(:) V(:) kz(:)]'; 
   if ~isreal(kk), warning('k is not real vector'); keyboard; end
   
   % the actual aoi for each (u,v) pair [deg]
   % list the aoi for each plane wave
   psi= acos(naoi' * kk); 
   
   % interpolate thin film response onto unique aoi's
   [psi_table, ipsi, itab] = unique(psi);
   Ttable = interp1(aoilist, RD(nl,:), psi_table, 'spline', 0);
   T = Ttable(itab);
      
   %T  = interp1(aoilist,RD(nl,:),psi,'spline',0); % if psi>aoi, set T = 0
   
   % apply response to angular spectrum
   T  = reshape(T,size(U)) .* S; 
  
   % apply defocus
   if abs(R_defocus) > eps,
       T = T .* exp(j*2*pi*(U.^2 + V.^2).*0.5.*R_defocus);
   end
   
   % propagation
   if abs(Z_propagate) > eps,
       T = T .* exp(j*2*pi*kz*Z_propagate);
   end
   
   % coupling efficiency
   % not allowing tilt to save the time of the fft
   %    NUM = abs(fft2(T.*S)).^2;
   %    TT0 = NUM(1,1);
   %    TTm = max2d(NUM);
   TStmp = T.*S;
   TT0 = abs(sum(TStmp(:))).^2; % = dc value of fft2(T.*S)
   eta0(nl) = TT0./SS0;
   %etam(nl) = TTm./SS0;
   etam(nl) = eta0(nl);
   
end
   
eta0 = eta0(:);
etam = etam(:);

    function [pol, aoi, Naoi, R_defocus, Z_propagate] = CheckOptions
        % default values
        pol = 'P';
        aoi = 0;
        Naoi = 400;
        R_defocus = 0;
        
        if exist('options','var'),
            if isfield(options,'Naoi'), Naoi = options.Naoi; end
            if isfield(options,'pol'), pol = options.pol; end
            if isfield(options,'aoi'), aoi = options.aoi; end
            if isfield(options,'R_defocus'), R_defocus = options.R_defocus; end
            if isfield(options,'Z_propagate'), Z_propagate = options.Z_propagate; end
        end

    end % CheckOptoins

end % main