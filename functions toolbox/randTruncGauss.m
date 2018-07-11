function x = randTruncGauss(N, mu, sig, w)
% x = randTruncGauss(N, mu, sig, w)
%
% truncated normal distribution
% N = size of array to return
% mu = mean of distribution
% sig = standard deviation of normal distribution
% w = truncated width factor, i.e. distribution is truncated at mu +/- w*sig

x = sig*randn(N);

while any(abs(x(:)) > w*sig),
    irep = abs(x) > w*sig;
    x(irep) = sig*randn(size(x(irep)));
end

x = x + mu;

return