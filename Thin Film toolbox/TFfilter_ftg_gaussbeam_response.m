function [eta0, etam, lamlist] =...
   TFfilter_ftg_gaussbeam_response(aoi0,itulam,deslam,beamdiam,pol)
% [eta0, etam, lamlist] = ...
%       TFfilter_gaussbeam_response(aoi0,itulam,deslam,beamdiam,pol)
%
%   aoi0 = nominal aoi to filter [deg]
%   itulam = nominal CWL at this aoi [um]
%   deslam = design wavelength for filter prescription [um]
%   beamdiam = beam waist diameter [um]
%   pol = 'P' or 'S'

constants;

% [um] beam waist radius
wx = beamdiam/2;  wy = beamdiam/2; 

% wavelength range for output
lam1 = itulam - 0.0015; lam2 = itulam + 0.0015;
Nlam = 211;

% AOI range for calculating filter response
Naoi = 400;
un = 6*itulam/(sqrt(2)*pi*wx); % max plane wave for calc filter
undeg = asin(un)/P; % convert to aoi
aoilist = linspace(aoi0-undeg,aoi0+undeg,Naoi); % aoi
if ~isreal(aoilist), warning('error aoilist'); keyboard; end

% for each aoi, get the response
TD = []; LM = []; , cwl = []; bw05 = []; bwadj = [];
for aoi = aoilist,
	
	[rd, td, lam] = ftg_spectrum(lam1,lam2,Nlam,...
		'aoi',num2str(aoi),'designwave',num2str(1000*deslam),'pol',pol);
	
	TD = [TD td(:)];  % lam rows X aoi columns
      
end

%figure, imagesc(aoilist,lam,abs(TD))

% for each wavelength, use angular spectrum to calculate coupling coeff.
lamlist = lam;
% filter normal
naoi = [sin(aoi0*P); 0; cos(aoi0*P)];

% create grid in angular spectrum domain
Nsig = 24; % umax: S(umax) = exp(-(Nsig^2)/2)
Nu   = 256;                   Nv   = 256;
umax = Nsig/(sqrt(2)*pi*wx);  vmax = Nsig/(sqrt(2)*pi*wy);
du   = 2*umax/Nu;             dv   = 2*vmax/Nv;
u    = [-Nu/2:Nu/2-1]'*du;    v    = [-Nv/2:Nv/2-1]'*dv;
[U, V] = meshgrid(u,v);

% calculate plane wave spectrum of gaussian beam
S = pi*sqrt(wx*wy)*exp(-(pi.^2)*((wx.*U).^2 + (wy.*V).^2));
SS0 = abs(sum(sum(S.^2))).^2;

for nl = 1:length(lamlist),
   lam  = lamlist(nl);

   % estimate filter response for each plane wave
   % propagation vector, u = sin(theta)/lam:
   kk = lam*[U(:) V(:) sqrt((1/lam).^2-U(:).^2-V(:).^2)]'; 
   if ~isreal(kk), warning('k is not real vector'); keyboard; end
   
   psi= acos(naoi' * kk)./P; % the actual aoi for each (u,v) pair [deg]
   T  = interp1(aoilist,TD(nl,:),psi,'spline',0); % if psi>aoi, set T = 0
   
   % apply response to angular spectrum
   T  = reshape(T,size(U)) .* S; 
  
   % coupling efficiency
   NUM = abs(fft2(T.*S)).^2;
   TT0 = NUM(1,1);
   TTm = max2d(NUM);
   eta0(nl) = TT0./SS0;
   etam(nl) = TTm./SS0;
   
end
   
eta0 = eta0(:);
etam = etam(:);

return