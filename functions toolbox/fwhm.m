function [width, rmax, xmax, errmsg] = fwhm(x,response)
% [width, rmax, xm] = fwhm(x,response)
% 
% width = fwhm
% rmax = max response, half-max is half of this.
% xmax   = interpolated sub-pixel location of the peak

% initialize return values to error
width = NaN;
rmax  = NaN;
xmax  = NaN;
errmsg = '';

if length(x) ~= length(response) || length(x) < 3,
    error('input vectors must have three points');
end

% use interpolated peak value
[x, response] = filterdata(~isnan(response) & ~isnan(x),x,response);
[rmaxtmp, xmaxtmp] = findpeakn(x,response,1);
if isempty(rmaxtmp),
    errmsg = 'no peak found';
    return
else
    [rmax, xmax] = deal(rmaxtmp, xmaxtmp);
end

% % correct for base line offset
% % average response excluding peak
% find valley to the left of the peak
[rtmp, xtmp] = filterdata(x < xmax, response, x);
% check that peak is not too close to end
if length(rtmp) <= 2, errmsg = 'peak too close to left end'; return, end 
ii = length(rtmp);
while ii > 1 && rtmp(ii-1) < rtmp(ii), ii = ii - 1; end
xleft = xtmp(ii);
% find valley to the right of the peak
[rtmp, xtmp] = filterdata(x > xmax, response, x);
if length(rtmp) <= 2, errmsg = 'peak too close to right end'; return, end %
ii = 1;
while ii < length(rtmp) && rtmp(ii+1) < rtmp(ii), ii = ii + 1; end
xright = xtmp(ii);

% average response excluding peak
rtmp = filterdata(x <= xleft | x >= xright, response);
response = response - mean(rtmp);
rmax = rmax - mean(rtmp);

% now set threshold to half max
[xf, rf] = filterdata( response >= 0.5*rmax, x, response);
[xedgelow, ii] = min(xf); redgelow = rf(ii);
[xedgehi,  ii] = max(xf); redgehi  = rf(ii);

% test if edge is at limits
if xedgelow <= min(x)
    %warning('lower limit');
    errmsg = [errmsg '; lower limit'];
    x1 = xedgelow;

% check if exact
elseif redgelow == 0.5*rmax,
    x1 = xedgelow;

else
    % interpolate
    [xflow, rflow] = filterdata( x < xedgelow, x, response );
    [xm, ii] = max(xflow); rm = rflow(ii);
    x1 = interp1([rm redgelow], [xm xedgelow], 0.5*rmax);
end

if xedgehi >= max(x),
    %warning('upper limit');
    errmsg = [errmsg '; upper limit'];
    x2 = xedgehi;
    
elseif redgehi == 0.5*rmax,
    x2 = xedgehi;
    
else,
    [xfhi,  rfhi]  = filterdata( x > xedgehi,  x, response );
    [xm, ii] = min(xfhi);  rm = rfhi(ii);
    x2 = interp1( [redgehi rm], [xedgehi xm], 0.5*rmax);
end

% ig = find(response > 0.5*rmax);
% if ig(1) > 1,
%     ii = [ig(1)-1 ig(1)];
%     x1 = interp1(response(ii), x(ii), 0.5*rmax);
% else,
%     x1 = x(ig(1));
% end
% if ig(end) < length(x),
%     ii = [ig(end) ig(end)+1];
%     x2 = interp1(response(ii), x(ii), 0.5*rmax);
% else,
%     x2 = x(ig(end));
% end

width = x2 - x1;

if isnan(width),
    keyboard;
end
