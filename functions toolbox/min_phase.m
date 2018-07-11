function [minp, gd, cd] = min_phase(lam,H,varargin)
% [phi, gd, cd] = min_phase(lam,H,options)
% lam = wavelength (must be sampled on a regular grid)
% H = amplitude response
% options: 'amplitude'  H is amplitude response (default)
%          'intensity'  H is intensity = amplitude^2
%          'dB'         H is 20*log10(amplitude)
%          'um'         lam is [um] (default)
%          'nm'         lam is [nm]
%          'THz'        lam is frequency [THz]
%          'Nfft'       followed by length of fft as string 
%                       (default = smallest power of 2 > length(lam)) 
%          'Ngrad'      order of gradient filter (default = 2)
%          'Fgrad'      cutoff freq of grad filter (default = 0.7 normalized)
%
% output: phi = [rad]
%         gd  = group delay [ps]
%         cd  = chromatic dispersion [ps/nm]

constants;
% set defaults
CC  = 1e6*C;     % convert to [nm . GHz]
lam = 1000*lam;  % convert from [um] to [nm] or [THz] to [GHz]
il  = abs(H);    % default value is amplitude
Nfft = 2^(ceil(log2(length(lam))));
resample_flag = 1; % default value

if ~isempty(varargin),
   if strmatch('amplitude',varargin),
      il = abs(H);
   elseif strmatch('intensity',varargin),
      il = sqrt(abs(H));
   elseif strmatch('dB',varargin),
      il = 10.^(H/20);
   else,
      il = abs(H);
   end
   
   if strmatch('nm',varargin),
      lam = 0.001*lam;
   end
   if strmatch('THz',varargin),
      ff = lam; % [GHz]
      dff = ff(2)-ff(1); % [GHz]
      lam = linspace(CC./ff(end),CC./ff(1),length(ff)); % [nm]
      resample_flag = 0;
   end
   
   infft = strmatch('Nfft',varargin);
   if infft, 
      if varargin{infft+1} > length(lam), 
         Nfft = varargin{infft+1};
      else,
         warning(['FFT length ' num2str(varargin{infft+1})...
               ' too small, using Nfft = ' num2str(Nfft)]);
      end
   end
   
   ingrad = strmatch('Ngrad',varargin);
   if ingrad, ngrad = sscanf(varargin{ingrad+1},'%d');
   else, ngrad = 2;
   end
   
   ifgrad = strmatch('Fgrad',varargin);
   if ifgrad, fgrad = sscanf(varargin{ifgrad+1},'%f');
   else, fgrad = 0.7;
   end
   
end

dlam = diff(lam);
dll = mean(diff(lam)); %lam(2)-lam(1);
ii = find(dlam<=0);
if ii, 
   warning('overlapping wavelengths');
   if ii == 1, lam(1) = lam(1)-dll;
   else, disp('don''t know what to do'); keyboard;
   end
end

% resample il onto frequency grid
if resample_flag,
   ff  = linspace(CC./lam(end),CC./lam(1),length(lam)); % [GHz]
   dff = ff(2)-ff(1);
   ilf = pchip(CC./lam,il,ff)';
else,
   ilf = il;
end
% call hilbert transform to get minimum phase
minp  = minphase(ilf,Nfft);

% group delay is gradient of phase (dff is [GHz])
%gdpf  = 1e3*gradient(minp,dff)/(2*pi); % [ps]
gdpf  = 1e3*gradfilt(minp,ngrad,fgrad)/(2*pi*dff); % [ps]

% resample GD back onto wavelength grid
gd = pchip(CC./ff,gdpf,lam);

% gradient in wavelength for chromatic dispersion (dll is [nm])
cd = gradient(gd,dll);

return

function minp = minphase(ilf,Nfft)
p0  = log(ilf);
p1  = ifft(p0-min(p0),Nfft);
p2  = [0 -sign([-Nfft/2+1:Nfft/2-1])]' .* p1;
p3  = real(j.*fft(p2));
minp = p3(1:length(ilf));
return

function df = gradfilt(f,n,fc)
b = firls(n,[0 fc fc+0.05 1],[0 1 0 0],'differentiator');
NN = floor(length(b)/2);
bnorm = -2*sum([-NN:-1].*b(1:NN));
df = filter(b,1,f)./bnorm;
df = df([floor(n/2)+1:end 1:floor(n/2)]);
return