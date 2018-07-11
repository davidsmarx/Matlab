function val = intgrt(x,y)
%val = intgrt(fun,xf,wf,x,y)
%given a function fun
%		 abcissas x
%      weights  wf
%output the approximate integral y = sum(wf(n)*fun(x(n)))

%get weights and abscissas for integration
global xf wf

if isempty([xf, wf])
   [xf, wf] = gausswts(-0.9,0.9,20);
   'initialized weights'
end

%pl = exp(i*k*(x*cos(t) + y*sin(t)));
[xx,yy] = meshgrid(x,y);
val = xx;
val(:)=wf'*exp(i*2*pi*(xf*(xx(:)') + sqrt(1 - xf.^2)*(yy(:)')));
