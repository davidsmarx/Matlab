n1 = 1;
n2 = 1.45;
nf = 1;
n3 = 1.45;
dlen = 1;    % [um] grating period
blen = 0.5;  % duty cycle
alen = 1;    % [um] grating depth
wavl = 1.55; % [um] wavelength
aoi  = 5*P;  % [rad]
gphi = -pi/2;% [rad] phase of grating
defo = 0;    % [waves] defocus aberation of incident beam
Nord = 10;   % # of diffraction orders in calculation

k0 = 2*pi./wavl;
k1 = n1.*k0;
k3 = n3.*k0;
ep2 = n2.^2;
epf = nf.^2;
Kg  = 2*pi./dlen;

% compute grating fourier coefficients
ii = [0:Nord-1];
alprow = blen.*(ep2-epf).*exp(-j*2*pi*ii.*(0.5*blen+gphi)).*sinc(pi*blen*ii);
alprow(1) = blen*ep2 + (1-blen)*epf;
alp = toeplitz(alprow);

% list of diffraction orders for calculation
iord = [floor(-(Nord-1)/2):floor((Nord-1)/2)];
i0   = floor(Nord/2)+1; % index of zero order
if length(iord)~=Nord, error('check #orders'); end

% compute relevant k's
kx  = k1.*sin(aoi) - iord.*Kg;
kz1 = conj(sqrt(k1.^2 - kx.^2));
kz3 = conj(sqrt(k3.^2 - kx.^2));
k2z0 = sqrt(k0.^2.*alprow(1) - (k1.*sin(aoi)).^2);

% compute matrix A for state equations
A00 = j*k2z0*eye(Nord);
A01 = j*eye(Nord);
A10 = j*k0.^2.*alp - j*diag(kx.^2);
A = [A00 A01; A10 A00];
clear A00 A01 A10;

[W, Mu] = eig(A);

% compute boundary condition matrix
mu = diag(Mu).';
M12 = [-eye(Nord) zeros(Nord); -diag(kz1) zeros(Nord)];
M21 = W .* (ones(2*Nord,1)*exp(alen.*(mu - j*k2z0)));
M22 = [zeros(Nord) -eye(Nord); zeros(Nord)  diag(kz3)];
M = [W M12; M21 M22];
clear M12 M21 M22;

% compute right hand side
rhs = zeros(4*Nord,1);
rhs(i0) = exp(-j*defo.*real(kz1(i0)));
rhs(Nord+i0) = -kz1(i0).*rhs(i0);

x = M \ rhs;

r = x(2*Nord+1:3*Nord);
t = x(3*Nord+1:4*Nord);

R = abs(r).^2.*real(kz1(:))./real(kz1(i0))
T = abs(t).^2.*real(kz3(:))./real(kz1(i0))

sum(R) + sum(T)