function s = gausspuls(t,T)
% s = gausspuls(t,T)
% return a gaussian shaped pulse with rms bandwidth = B = 1/T, T = sqrt(2)*sigma, 
% sigma = time domain variance
% t = time is a real 1-d vector

sigma = T / (2*sqrt(log(10)));
sigma2 = sigma.^2;
s = exp( -(t.^2)/(2*sigma2) );

return