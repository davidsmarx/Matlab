function [Pzz, fx, fy] = psd2d(Z, nfft, Fsx, Fsy, window)
%[Pzz, fx, fy] = psd2d(Z, nfft, Fsx, Fsy, window)

if nargin == 0 | nargin > 5, 
   error('usage: [Pzz, fx, fy] = psd2d(Z, nfft, Fsx, Fsy, window)\n');
elseif nargin <= 4,
   if nargin <= 3,
      Fsy = 1;
      if nargin <= 2,
         Fsx = 1;
         if nargin == 1,
            nfft = 128;
         end
      end
   end
   window = hanning(nfft);
end

   
[Ny Nx] = size(Z);
Nw = length(window);
if Nw > nfft, error('window is longer that FFT length'); end

w = zeros(nfft,1);
w(1:Nw) = window(1:Nw);
W = w * w';
disp(['size of window matrix = ' num2str(size(W))]);

Pzz = zeros(nfft);

cnty = 0;
for ny = 1:floor(nfft/2):Ny-nfft,
   cntx =  0;
   for nx = 1:floor(nfft/2):Nx-nfft,
      
      ztmp = Z(ny:ny+nfft-1,nx:nx+nfft-1) .* W;
      pp = fft2( ztmp - mean(mean(ztmp)) );
      Pzz = Pzz + pp.*conj(pp);
      cntx = cntx + 1;
   end
   cnty = cnty + 1;
end
disp(['did ' num2str(cntx) ' x ' num2str(cnty) ' 2d fft']);

Pzz = fftshift(Pzz) ./ (cnty*cntx);

fx = Fsx.*[ -0.5 : 1/nfft : 0.5 - 1/nfft ];
fy = Fsy.*[ -0.5 : 1/nfft : 0.5 - 1/nfft ];

return
