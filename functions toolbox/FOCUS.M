function val = focus(x,y,na)
%val = focus(x,y,na)
%integrates exp(i*2*pi*(u*x + sqrt(1-u^2)*y)) du from -na to na
%using gaussian quadrature

%get weights and abscissas for integration
global xf wf g_na

if or( isempty([xf, wf]) , g_na~=na )
   [xf, wf] = gausswts(-na,na,20);
   g_na = na;
   'initialized weights'
end

%pl = exp(i*k*(x*cos(t) + y*sin(t)));
[xx,yy] = meshgrid(x,y);
val = xx;
val(:)=wf'*exp(i*2*pi*(xf*(xx(:)') + sqrt(1 - xf.^2)*(yy(:)')));
