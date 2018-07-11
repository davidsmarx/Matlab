function [resp, xr] = calc_response(cdata,x,y)
%[resp, xr] = calc_response(cdata,x,y)

[cm nr nc] = max2d(cdata);
nplot = ([nc-50:nc+50]);

resp = decibel(cdata(nr,nplot));
xr = x(nplot);

%figure, plot(xr,resp,'-o'), grid

return