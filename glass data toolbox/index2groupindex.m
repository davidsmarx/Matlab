function [group_index, Sout] = index2groupindex(Sin, lami)
% [group_index, Sout] = index2groupindex(Sin, lami)
%
% Sin.wavelength
% Sin.index
%
% lami = wavelength grid for calculating output

unitsdefinitions;
constants;

if isfield(Sin,'wavelength'),
    % optical frequency from the table:
    f_n = C./Sin.wavelength;
    f_nmin = min(f_n);
    f_nmax = max(f_n);
    lam_n = Sin.wavelength;
    lammin = min(Sin.wavelength);
    lammax = max(Sin.wavelength);
else,
    error('usage');
end

% index of refraction from the table:
n = Sin.index;

% fit a polynomial to the index vs. optical frequency:
polyorder = min(4, floor(length(f_n)/2));
[p_n, s, mu] = polyfit(f_n, n, polyorder);

% plot index vs. wavelength for table and polynomial fit
lamplot = linspace(lammin, lammax)';
% figure, plot(lamplot/UM, polyval(p_n, C./lamplot,[],mu), '-', lam_n/UM, n, 'o'), grid
% xlabel('Wavelength (\mum)'), ylabel('Index')

freqplot = linspace(f_nmin,f_nmax)';
% figure, plot(freqplot/THZ, polyval(p_n, freqplot, [], mu), '-', f_n/THZ, n, 'o'), grid
% xlabel('Frequency (THz)'), ylabel('Index')

% group index - lambda
p_dndf = p_n(1:end-1).*(polyorder:-1:1)./mu(2);
funGroupindex = @(f)(polyval(p_n, f, [], mu) + f.*polyval(p_dndf, f, [], mu));
% figure, plot(freqplot/THZ, funGroupindex(freqplot)), grid,
% xlabel('Frequency (THZ)'), ylabel('Group Index')

% outputs:
group_index = funGroupindex(C./lami);
Sout.funGroupindex = funGroupindex;
Sout.polyfit = struct('p',p_n,'s',s,'mu',mu);
Sout.freqplot = freqplot;
Sout.groupindexplot = funGroupindex(freqplot);
