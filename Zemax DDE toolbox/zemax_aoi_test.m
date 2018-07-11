[Lens, zchan] = z_Getlens;

wedge = 9.05; % deg

% set wedge angle
status = z_setpar(zchan,Lens.WEDGE_ANGLE_1,3,wedge);
status = z_setpar(zchan,Lens.WASHER,2,tan(wedge*pi/180));

% SET CRYSTAL ORIENTATION
% par2 = x-cosine, par3 = y-cosine, par4 = z-cosine
% par1 = 0 => ordinary, par1 = 1 => extraordinary
status = z_setpar(zchan,Lens.WEDGE_1_IN,[2 3 4],[0 1 0]);
status = z_setpar(zchan,Lens.WEDGE_2_IN,[2 3 4],[1 0 0]); 

h = [0 0]; p = [0 0];
config_f_eo = 1;
config_f_oe = 2;

for tilt = -10:5,
   
   status = z_setpar(zchan,Lens.TILT_WHOLE_PRISM,3,tilt);
   
   for config = [config_f_eo config_f_oe],
      status = z_setconfig(zchan,config);
      
      [vc, raytest, raycostest, sn, rayI(config)] = ...
         z_gettrace(zchan,1,0,Lens.OUTPUT.nsurf,h,p);
      [vc, rayimage, raycosimage, sn, rayItmp] = ...
         z_gettrace(zchan,1,0,Lens.IMAGE.nsurf,h,p);
      y(config) = raytest(2);
      yc(config) = raycostest(2);
      yf(config) = rayimage(2);
      
   end

   thsplit = (90/pi)*(asin(yc(config_f_eo))-asin(yc(config_f_oe)));
   thbis = (90/pi)*(asin(yc(config_f_eo))+asin(yc(config_f_oe)));
   zsplit = (y(config_f_oe) - y(config_f_eo)) /...
      (tan(asin(yc(config_f_eo))) - tan(asin(yc(config_f_oe))));
   fibsep = abs(yf(config_f_oe)-yf(config_f_eo));
   
   fprintf('%f, %f, %f, %f, %f\n',tilt,thsplit,thbis,zsplit,fibsep);
  
end % tilt


rc = ddeterm(zchan);
