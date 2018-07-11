function y = rms(x)
% y = sqrt(mean(x.^2));
% 
% x is a vector, don't use with x a matrix
y = sqrt((x(:)'*x(:))./length(x(:)));
