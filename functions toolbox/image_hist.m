function [histogram, hist_rang] = image_hist(win1,N)
% [hist, hist_range] = image_hist(image,N)
% image is a real image
% N is the number of bins for the histogram [default = 100]
% or if N is a vector, it is the bins

if nargin == 0,
   error('usage: [hist, hist_range] = image_hist(win1,N)');
end

if nargin == 1,
   N = 100;
end

if length(N) == 1, % N is scalar
   maxs = max(max(win1));
   mins = min(min(win1));
   bins = (maxs - mins)/(N-1);
   hist_rang = (mins:bins:maxs)';
   
else
   hist_rang = N;
end

histogram = sum(hist(win1,hist_rang),2);
