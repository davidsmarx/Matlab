function y = rms(x, idim)
% y = rms(x)
% y = rms(x, idim)
%
% y = sqrt(mean(x.^2, idim));

%y = sqrt((x(:)'*x(:))./length(x(:)));

if exist('idim','var'),
    y = sqrt(mean(x.^2, idim));
else
    y = sqrt(mean(x.^2));
end
