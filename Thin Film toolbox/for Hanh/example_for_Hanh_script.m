clear all;

% path for index of silicon
addpath('glass data toolbox');

% units:
U = CConstants;

% create a reflected spectrum from a Si wafer with 1um oxide on the front
% side

% thickness of materials
SiThick = 775*U.UM;
OxThick = 1.0*U.UM;

% index of materials
[lamSi, nSi] = IndexOfSi;
nOx = 1.46; % can replace with a table of wavelength and index if you want to include dispersion

% WTS wavelengths:
Nlam = 2048;
lam = linspace(1480,1620, Nlam)'*U.NM;

% calculate the reflected intensity at each wavelength => spectrum
R = zeros(Nlam,1);
for ii = 1:Nlam,
    % interpolate the index of Si at this wavelength
    nSi_ii = interp1(lamSi, real(nSi), lam(ii));
    R(ii) = thin_film_filter_2([1 nOx nSi_ii 1], [OxThick SiThick], 0, lam(ii));
end

% calculate thickness from the spectrum
[thick, amp, S] = spectrum2thickness(lam, R, [lamSi real(nSi)]);

% display the results
fprintf('calculated thickness = %.3fum\n', thick/U.UM);

% plot the spectrum
figure, plot(lam/U.NM, R), grid
xlabel('Wavelength (nm)'), ylabel('R')

% plot the A-Scan
figure, plot(S.zd/U.UM, S.mag_td), grid
xlabel('Thickness (\mum)')
ylabel('Amplitude')
