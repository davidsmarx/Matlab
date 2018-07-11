%try_sphere_fit
%
% IB 
%
% revival of a 13 years old code


  % Create data for a circle + noise
  
  ph = linspace(0,2*pi,20)';
  th = linspace(0,pi,20)';
  [PH, TH] = meshgrid(ph,th);
  R=1.1111111;
  sigma = R/10;
  z = R.*cos(TH) + randn(size(TH))*sigma;
  y = R.*sin(TH).*sin(PH) + randn(size(TH))*sigma;
  x = R.*sin(TH).*cos(PH) + randn(size(TH))*sigma;
  
  plot3(x,y,z,'o'), title(' measured points')
  pause(1)
   
  % reconstruct circle from data
  [xc,yc,zc,Re,a] = spherefit(x,y,z);
  xe = Re*sin(TH).*cos(PH)+xc;
  ye = Re*sin(TH).*sin(PH)+yc;
  ze = Re*cos(TH)+zc;
  
  plot3(x(:),y(:),z(:),'o',[xe(:);xe(1)],[ye(:);ye(1)],[ze(:);ze(1)],'-.',...
      R*sin(TH(:)).*cos(PH(:)),R*sin(TH(:)).*sin(PH(:)),R*cos(TH(:)),'-'),
  title(' measured fitted and true circles')
  legend('measured','fitted','true')
  fprintf('center (%g , %g, %g );  R=%g\n',xc,yc,zc,Re)
  xlabel x, ylabel y
  axis equal
     

   
  


