function [thickness, amplitude, spectrum] = spectrum2thickness(varargin)
% [thickness, amplitude, spectrum] = spectrum2thickness(lam, level, nfilm)
% [thickness, amplitude, spectrum] = spectrum2thickness(spectrum, nfilm)
%
% spectrum is a struct:
%     frequency or wavelength
%     level
%     minthickness (optional)
%     maxthickness (optional)
%     minwavelength (optional)
%     maxwavelength (optional)
%
% nfilm = scalar group index of refraction, or
% nfilm = [wavelength index] array, wavelength range must cover data
%
% outputs:
%   spectrum.fr
%   spectrum.mag_fr
%   spectrum.td
%   spectrum.zd
%   spectrum.mag_td

C = 2.99792458e8; % [vacuum m/s]

if nargin == 3,
    % old method:
    [spectrum.wavelength, spectrum.level, nfilm] = deal(varargin{:});
else
    [spectrum, nfilm] = deal(varargin{:});
end

% if input is wavelength, resample spectrum onto frequency grid
if isfield(spectrum,'wavelength')
    
    % check for optional min/max wavelength
    if isfield(spectrum,'minwavelength'),
        [spectrum.wavelength, spectrum.level] = filterdata(...
            spectrum.wavelength > spectrum.minwavelength, ...
            spectrum.wavelength, spectrum.level);
    end
    if isfield(spectrum,'maxwavelength'),
        [spectrum.wavelength, spectrum.level] = filterdata(...
            spectrum.wavelength < spectrum.maxwavelength, ...
            spectrum.wavelength, spectrum.level);
    end
    
    N = 2048; % choose length of resampled spectrum
    ff = C./spectrum.wavelength;
    fr = linspace(min(ff),max(ff),2048)';

    sr = interp1(ff, spectrum.level, fr);
    % figure, plot(fr/THZ, sr), grid,...
    %     xlabel('Frequency (THz)'), ylabel('Intensity'),...
    %     title(pwd2titlestr(filename))
    
elseif isfield(spectrum,'frequency')
    fr = spectrum.frequency;
    sr = spectrum.level;
    N = length(fr);
    
    % check if fr is ascending or descending
    if fr(1) > fr(end),
        fr = flipud(fr);
        sr = flipud(sr);
    end
    % check lengths
    if ~isequal(length(fr), length(sr)),
        error('input arrays must have same length');
    end
else
    error('usage');
end

% check spectrum sr for NaN's because one NaN ruins the whole fft
if any(isnan(sr)),
    error('input array has a NaN');
end

% if nfilm is an array, do the frequency scaling method
if ~isscalar(nfilm)
    % do checks
    
    % interpolate index of refraction onto data frequency grid
    disp_lam = nfilm(:,1);
    disp_n   = nfilm(:,2);
    disp_freq = C./disp_lam;
    disp_indx_r = interp1(disp_freq, disp_n, fr);
    
    % multiply index * frequency
    fr_n = fr .* disp_indx_r;
    
    %     figure, plot(fr_n/1e12, sr), grid
    %     xlabel('Frequency * n'), ylabel('Spectrum Amplitude')
    
    % resample spectrum onto even (index*frequency) grid
    fr_n_r = linspace(fr_n(1),fr_n(end),length(fr_n))';
    sr_r   = interp1(fr_n, sr, fr_n_r);
        
    % all done, copy resampled spectrum to working spectrum, and set
    % nfilm = 1 since we don't need to correct anymore
    fr = fr_n_r;
    sr = sr_r;
    nfilm = 1;
    
end
% remove low order curvature, due to system's spectral response
[p, s, mu] = polyfit(fr, sr, 3);
sr = sr - polyval(p, fr, [], mu);

% apply window
sr = sr.*hanning(N);

% 
df = mean(diff(fr));

% transform to time delay domain, pad fft by a lot for better peak finding
Nfft = 11*N;
dt = 1./(Nfft*df);
ftsp = fft(sr, Nfft);
%figure, plot(dt*(0:N/2-1)',abs(ftsp(1:N/2)),'-o'), grid

tt = dt*(0:floor(Nfft/2)-1)';
ftsp = ftsp(1:floor(Nfft/2));

% determine thickness from peak
if isfield(spectrum, 'minthickness'), mintt = spectrum.minthickness * (2 * nfilm ./ C);
else mintt = tt(10); end
if isfield(spectrum, 'maxthickness'), maxtt = spectrum.maxthickness * (2 * nfilm ./ C);
else maxtt = tt(end); end

[ttuse, ftspuse] = filterdata(tt >= mintt & tt <= maxtt, tt, ftsp);
[fm, tpeak] = findpeakfftinterp(ttuse, abs(ftspuse), 1);

thickness = C .* tpeak ./ (2 * nfilm);
amplitude = fm;

spectrum.fr = fr;
spectrum.mag_fr = sr;
spectrum.td = tt';
spectrum.zd = C .* tt ./ (2 * nfilm);
spectrum.mag_td = abs(ftsp);
spectrum.complex_td = ftsp;

%fprintf('thickness = %fum, tpeak = %fps\n',thickness/UM,tpeak/PS);
