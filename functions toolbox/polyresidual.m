function [yres, p, s, mu] = polyresidual(x,y,n)
% [yresidual, p, s, mu] = polyresidual(x,y,n)

[xtmp, ytmp] = filterdata(~isnan(x) & ~isnan(y), x, y);

[p, s, mu] = polyfit(xtmp, ytmp, n);

yres = y - polyval(p, x, [], mu);

