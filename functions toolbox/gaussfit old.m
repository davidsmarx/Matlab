function [a, sig] = gaussfit(y,x)
% [a, sig] = gaussfit(y,x)
%
% to place center of mass at the origin:
% com = inline('x - mean(x.*u)./mean(u)','x','u');



N = length(y);
if nargin == 0, error('usage: [a, sig] = gaussfit(y,x)'); end
if nargin == 1, x = [-N/2:N/2-1]'; end
if nargin >= 2, x = reshape(x,size(y)); end

sig0 = sqrt(sum(x.^2.*y)./sum(y));

sig = fminsearch(@fitMF,sig0,[],x,y);
a   = calcamp(x,y,sig);

return


function a = calcamp(x,y,sig)

gtmp = exp(-(x./sig).^2);

a = sum(y.*gtmp)./sum(gtmp.^2);

return


function I = fitMF(sig,x,y)

a = calcamp(x,y,sig);
I = sum((y - a.*exp(-(x./sig).^2)).^2);

return