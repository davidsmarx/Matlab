function [fm, xm] = findpeakfftinterp(x,y,n)
% [fm, xm] = findpeakfftinterp(x,y)
% [fm, xm] = findpeakfftinterp(x,y,n)
% [fm, xm] = findpeakfftinterp(x,y,n,padfact)
%
% assumes x is on a evenly spaced grid
% n (default = 1) -- find the n highest peaks
%    [fm(1:n), xm(1:n)] are in descending order of peak amplitude
% padfact (default = 5) -- (padfact - 1) points are interpolated between
%    each input sample

if ~exist('n','var') || isempty(n),
    n = 1;
end

N = length(y);

if ~exist('padfact','var') || isempty(padfact),
    padfact = 5;
end

Y = fft(y(:));
Y = fftshift(Y);
Yp = [zeros(padfact*N+1,1); Y; zeros(padfact*N-1,1)];

yp = abs((2*padfact+1)*ifft(fftshift(Yp)));

dx = mean(diff(x));
dxp = dx./(2*padfact + 1);
xp = x(1) + dxp*(0:N*(2*padfact+1)-1)';

if n > 1,
    
    % re-written using Matlab findpeaks()
    [pks, locs] = findpeaks(yp, 'sortstr', 'descend');
    [fm, xm] = deal(zeros(n,1));
    for ii = 1:min(n, length(pks)),
        [fm(ii), xm(ii)] = ...
            findpeakn(xp(locs(ii) + [-1 0 1]), yp(locs(ii) + [-1 0 1]), 1);
    end
    
    %%%% old method
    %     [ysort, isort] = sort(yp,1,'descend');
    %
    %     % go one-by-one through the highest psd values, and select those that are
    %     % peaks. vpsd(ii) is a peak if vpsd(ii) > vpsd(ii-1) && vpsd(ii) >
    %     % vpsd(ii+1)
    %     nfound = 0;
    %     ip = 0;
    %     [fm, xm] = deal(zeros(n,1));
    %     while nfound < n,
    %         ip = ip+1;
    %         itest = isort(ip);
    %         if itest > 1 && itest < length(yp),
    %             if yp(itest) > yp(itest-1) && yp(itest) > yp(itest + 1),
    %                 nfound = nfound + 1;
    %                 %                 ypeak(nfound) = yp(itest);
    %                 %                 xpeak(nfound) = xp(itest);
    %                 [fm(nfound), xm(nfound)] =...
    %                     findpeakn(xp(itest + [-1 0 1]), yp(itest + [-1 0 1]), 1);
    %             end
    %         end % if not an endpoint
    %     end
 
    
else,
    [fm, xm] = findpeakn(xp, yp, 1);
end

