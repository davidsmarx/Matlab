function [mtf, ff] = calcMTF(Iroi,options)
% [mtf, ff] = calcMTF(Iroi,options)
%
% calculate the MTF from the image of an edge (edge response)
%
% Iroi = grayscale image, or filename of an image, or [] for interactive
%
% fields for options struct:
%    method = 'spatial' (default), 'frequency'
%    direction = 'horizontal' (default), 'vertical'
%    bgradfilt = FIR filter coefficients for the gradient filter used in
%         'spatial' method. (default = [0.5 -0.5])
%    dx = pixel spacing (default = 1)

if nargin==0 | isempty(Iroi),
    [filename, pathname, filterindex] = uigetfile({'*.jpg;*.tif'},'Select Image File');
    I = imread([pathname filename]);
    [Iroi, x, y, h_fig] = getroi(I);
    
elseif isstr(Iroi),
    Iroi = imread(Iroi);
end
    
if ~exist('options','var'), options = struct; end
options = SetOptions(options);

switch lower(options.direction)
    case 'horizontal'
        Iroi = Iroi';
    case 'vertical'

    otherwise
        error('unrecognized direction option');
end

switch lower(options.method)
    case 'spatial'
        [mtf, ff] = spatialMTF(Iroi,options);
    case 'frequency'
        [mtf, ff] = freqMTF(Iroi,options);
    otherwise
        error('unrecognized method option');
end

return % calcMTF
        
function optout = SetOptions(optin)

optout = optin;
if ~isfield(optout,'method'), optout.method = 'spatial'; end
if ~isfield(optout,'direction'), optout.direction = 'horizontal'; end
if ~isfield(optout,'bgradfilt'), optout.bgradfilt = [1 -1]; end
if ~isfield(optout,'dx'), optout.dx = 1; end

return % SetOptions

function [mtf, ff] = spatialMTF(Iroi,options)

% get power of 2
[Ny, Nx] = size(Iroi);
Nfft = 2.^ceil(log2(Ny));
% calculate gradient (edge response to line response)
droi = filter(options.bgradfilt,1,double(Iroi));
droi(1:length(options.bgradfilt)-1,:) = 0; % (lack of) initial conditions
% mtf = abs(mean(fft(droi,Nfft),2));
mtf = mean(abs(fft(droi,Nfft)),2);

mtf = fftshift(mtf)./max(abs(mtf));
ff = [-Nfft/2:Nfft/2-1]'/Nfft/options.dx;

return % spatialMTF

function [mtf, ff] = freqMTF(Iroi,options)

% get power of 2
[Ny, Nx] = size(Iroi);
Nfft = 2.^ceil(log2(Ny));

I0 = mean(Iroi(end-9:end,1)) - mean(Iroi(1:10,1));

for ii = 1:Nx,
    atmp = double(Iroi(:,ii));
    ltmp = linspace(atmp(1),atmp(end),Ny)'; % subtract line so edges = zero
    mtmp(:,ii) = fftshift(fft(atmp-ltmp,Nfft));
end
ii = 2*pi*[-Nfft/2:Nfft/2-1]'./Nfft;
mtf = abs(ii.*mean(mtmp,2))./I0;
mtf(Nfft/2+1) = 1;

ff = [-Nfft/2:Nfft/2-1]'/Nfft/options.dx;


return % freqMTF



