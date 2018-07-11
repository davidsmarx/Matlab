function F = czt2d(f, dx, dy, dfx, dfy, inverse)
% F = czt2d(f, dx, dy, dfx, dfy, inverse)
%
% F = fftshift(czt2d(fftshift(f), 1, 1, 1/N, 1/N, false))
%   <=> 
% F = fft2(f)
%
% F = (1/N)*(1/N)*fftshift(czt2d(fftshift(f), 1, 1, 1/N, 1/N, true))
%  <=>
% F = ifft2(f)


if ~exist('inverse','var') | isempty(inverse),
    inverse = false;
end

[nyin, nxin] = size(f);
nyout = nyin; nxout = nxin; % force same size for now

F = zeros(nyout,nxout);

%! i.e.  Nfft = 2^ceil(log2(2*Ndata)) = 2^(1+ceil(log2(Ndata)))
nfftx = 2.^(1 + ceil(log2(nxin))); % Lcols
nffty = 2.^(1 + ceil(log2(nyin))); % Lrows

alphax = pi*dfx*dx;
alphay = pi*dfy*dy;

% transform each column
n1 = nyin/2;
k1 = nyout/2;

nymax = max( nyin, nyout );
nxmax = max( nxin, nxout );

ZA = LocalZPower( alphay, nyin, 2*k1, n1*k1 );
ZB = LocalZPower( alphay, nyout, 2*n1, n1*k1 );
Z_CHIRP = LocalZChirp ( alphay, nymax, nffty );

fft_Z_chirp = zeros(nffty,1);
fft_Z_chirp(1:nyout) = Z_CHIRP(1:nyout);
fft_Z_chirp(nffty-nyin+2:nffty) = flipud(Z_CHIRP(2:nyin));

fft_Z_chirp = fft(fft_Z_chirp);

for k = 1:nxin,
    F(:,k) = LocalCZT1D(f(:,k), nyin, nyout, nffty, ZA, ZB, fft_Z_chirp, inverse);
end

% transform each row
n1 = nxin/2;
k1 = nxout/2;

ZA = LocalZPower( alphax, nxin, 2*k1, n1*k1 );
ZB = LocalZPower( alphax, nxout, 2*n1, n1*k1 );
Z_CHIRP = LocalZChirp( alphax, nxmax, nfftx );

fft_Z_chirp = zeros(nfftx,1);

fft_Z_chirp(1:nxout) = Z_CHIRP(1:nxout);
fft_Z_chirp(nfftx-nxin+2:nfftx) = flipud(Z_CHIRP(2:nxin));

fft_Z_chirp = fft(fft_Z_chirp);

for jj = 1:nxout,
    F(jj,:) = LocalCZT1D(F(jj,:), nxin, nxout, nfftx, ZA, ZB, fft_Z_chirp, inverse);
end
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    function z_phasor = LocalZPower(alpha, n, two_k_half, nk1)
        k = (1:n)';
        z_phasor = exp( -j*alpha*( (k - two_k_half - 1).*(k-1) + nk1) );
    end

    function z_chirp = LocalZChirp(alpha, n, L)
        factor = 1./L;
        k = (1:n)';
        z_chirp = factor * exp( j.*alpha.*(k-1).^2 );
    end

    function czt = LocalCZT1D( array, numIn, numOut, L, za, zb, fft_z_chirp, invczt)

        array = array(:); % 
        u = zeros(L,1);
        czt = zeros(numOut,1);
        
        if invczt
            u(1:numIn) = conj(array) .* za;
        else
            u(1:numIn) = array .* za;
        end
        
        u = fft(u);
        u = u .* fft_z_chirp;
        % u = ifft(u);
        u = conj(fft(conj(u)));
        
        czt = u(1:numOut) .* zb;
        
        if invczt
            czt = conj(czt);
        end
    end % LocalCZT1D

end

