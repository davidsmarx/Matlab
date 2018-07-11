function [fm, xm, pp, rmsres] = findpeakn(xi,yi,n)
% [fm, xm, pp, rmsres] = findpeakn(xi,yi,n)
%
% quadratic fit to interpolate peak value
% xi = abscissa vector
% yi = vector of values at each abscissa
% n  = number of sample points around the maximum to use in quadratic fit.
% The points used are imax + [-n:n]; thus a total of 2*n+1 points are used.
%
% fm is the value of the peak
% xm is the location of the peak
% pp are the polynomial coefficients: pp(1)*x^2 + pp(2)*x + pp(3) = 0
% rmsres is the rms of the residual between the polynomial fit and the (xi,
% yi)

[ymax, imax] = max(yi);

% check that peak is not at an end point
if imax == 1,
    %warning('peak is at left end point');
    fm = yi(1);
    xm = xi(1);
    pp = [];
        
    return
elseif imax == length(yi),
    %warning('peak is at right end point');
    fm = yi(end);
    xm = xi(end);
    pp = [];
        
    return
end

ii = imax + [-n:n];
xx = xi(ii) - xi(imax);
[pp, S, mu] = polyfit(xx,yi(ii),2);

xt = -pp(2)/2/pp(1); xt = mu(1) + mu(2)*xt;
fm = polyval(pp,xt,[],mu);
xm = xt + xi(imax);

residual = yi(ii) - polyval(pp, xx, [], mu);
rmsres = rms(residual);

%%%% old way %%%%
% pp = polyfit(xx,yi(ii),2);
% 
% xt = -pp(2)/2/pp(1);
% fm = polyval(pp,xt);
% xm = xt + xi(imax);


return
