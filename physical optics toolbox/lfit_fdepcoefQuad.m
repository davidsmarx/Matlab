function [xx,residual,fval]=lfit_fdepcoefquad(xdata,ydata)
% [xx,residual,fval]=lfit_fdepcoefquad(xdata,ydata)
%ydata: column function to be linear and radial quadratic fitted
%xdata: fitting variables(positions) X:xdata(:,1) Y:xdata(:,2)
%xx: [DC XlinearSlope YlinearSlope radialQuadCoef]
%fval:fitted function
%resudual: residual after fitting
%by tjs 1-27-03
%update 3-25-04 by tjs

nlen=length(xdata);
Aones=ones(nlen,1);
A=[Aones xdata(:,1) xdata(:,2) xdata(:,1).^2+xdata(:,2).^2];


%Z = pinv(A)*data;%this choice is slower by 50%!!!!
P = inv(A'*A);
xx = P*A'*ydata;

  
  fval=A*xx;
	residual=ydata-fval;
   
   
   return
