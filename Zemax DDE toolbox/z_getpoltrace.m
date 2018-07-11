function [rayI, E] = z_getpoltrace(zchan,wave,mode,nsurf,h,p,Exy)
% [rayI, E] = z_getpoltrace(zchan,wave,mode,nsurf,h,p,Exy)
% wave = wavelength #
% mode: 0 = real, 1 = paraxial
% nsurf = surface to which ray is traced
% h = [hx, hy]
% p = [px, py]
% Exy = [Ex, Ey]

Exy = Exy(:);
Exy = Exy./sqrt(Exy'*Exy);

cmdstr = ['GetPolTrace,'...
      num2str(wave) ','...
      num2str(mode) ','...
      num2str(nsurf) ','...
      sprintf('%d,%d,%d,%d,',h,p)...
      sprintf('%e,%e,%e,%e',abs(Exy),angle(Exy)*180/pi)];

rays = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf])';
% rays(1) = error
% rays(2) = intensity
% rays(3:5) = real(Ex,Ey,Ez)
% rays(6:8) = imag(Ex,Ey,Ez)

if rays(1) > 0, 
   warning(['ray missed surface #' num2str(rays(1))]); keyboard;
elseif rays(1) < 0,
   warning(['TIR at surface #' num2str(-rays(1))]); keyboard;
end

rayI = rays(2);
E = rays(3:5) + j*rays(6:8);

return


