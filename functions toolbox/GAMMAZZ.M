function [f] = gammazz(z);
% GAMMAZZ
%        Complex Gamma function valid in the entire complex
%        plane. This routine uses the reflection formula
%        and the GAMMAZ routine to provide valid results
%        for all z. This routine should be used in place
%        of GAMMAZ because the latter is inaccurate for
%        real(z)<0. This routine uses the GAMMAZ routine.
%
%        For real arguments GAMMAZZ gives virtually the
%        same results as MATLAB's GAMMA routine.
%
%   algorithm: numerical recipes 6.1.5, "gammln"
%   Ref: "Special Functions ..." by L.C. Andrews
%
%   Tested under version 5.2.0
%
%   usage: [f] = gammazz(z)
%   z can be any size
%
%see also   GAMMA GAMMAZ GAMMAINC


%Paul Godfrey
%pjg@mlb.semi.harris.com


twopi=pi+pi;
[m n] = size(z);
z=z(:);
f = 0.*z; % reserve space in advance
r=f+1/2;  % reflect about 1/2



  zz =z;
  R=find(real(zz)>=r);
  L=find(real(zz)< r);


if length(R)>0;
  gp=gammaz(zz(R));
  ff(R)=gp;
end
  
if length(L)>0;
  zzz=zz(L);
  xx=pi.*real(zzz);yy=pi.*imag(zzz);
  xx=mod(xx+twopi,twopi);
  xp=mod(xx+pi,twopi);zp=xp+i*yy;
  xn=    xx-pi       ;zn=xn+i*yy;
% (sin(z+pi)+sin(z-pi))/(-2)=sin(z)
  sinz=(sin(zp)+sin(zn))./(-2);
  rz=real(zzz);iz=imag(zzz);
  s=find(rz==round(rz) & iz==0);
  sinz(s)=0.0;
  %gn = zeros(size(zzz));
  %gn(s) = Inf;
  %gn(~s) = pi./(sinz(~s).*gammaz(1.-zzz(~s)));
  gn = pi./(sinz.*gammaz(1.-zzz));
  ff(L)=gn;
end


k=find(~isfinite(ff));
ff(k)=Inf;


f=reshape(ff,m,n);


return


% A demo of this function is
x=-4:0.125:4;
y=x;
[X,Y]=meshgrid(x,y);
z=X+i*Y;
g=gammazz(z);
%Infinity doesn't plot well
%so for plotting purposes
%lets use 10
p=find(isinf(g));
g(p)=10;


mesh(abs(g));
rotate3d;
