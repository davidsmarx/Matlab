function [a, sig, x0, yp] = gaussfit(y,x)
% [a, sig, x0, yp] = gaussfit(y,x)
% 
% mse fit of log(y) to log(a*exp(-(x-x0)^2/(2*sig^2)))
%
% yp = fit gaussian evaluated on x
%
% for gaussian fit not using logarithms, see 'gaussfit old.m'

N = length(y);

if ~exist('x','var'), x = [1:N]'; end

% y must be strictly > 0
yoffset = -min(y) + 1e-9*range(y);
y = y + yoffset;

%%%%%%%%%%%%%%%% merit function for searching for x0
function [val, a, sig] = MF(x0)

    t = x - x0;

    % solve for amp and sig in log
    pp = [ones(N,1) -t.^2] \ log(y);
    a = exp(pp(1));
    sig = sqrt(0.5/pp(2));

    yp = a.*exp(-t.^2./(2*sig^2));
    yd = y - yp;
    val = yd'*yd;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% estimate offset
%x00 = mean(x.*y)./mean(y);
[ymax, imax] = max(y);
x00 = x(imax);

[x0, fval, exitflag, output] = fminsearch(@MF, x00);

[val, a, sig] = MF(x0);

yp = yp - yoffset;

end