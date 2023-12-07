function [vigcode, raypos, raycos, surfnorm, rayI] = ...
   z_gettracedirect(zchan, wave, mode, nsurf_start, nsurf_stop, raypos_start, raycos_start)
% [vigcode, raypos, raycos, surfnorm, rayI] = ...
%      z_gettracedirect(zchan, wave, mode, nsurf_start, nsurf_stop, raypos_start, raycos_start)
% wave = wavelength #
% mode: 0 = real, 1 = paraxial
% nsurf_start = nsurf for surface to launch ray
% nsurf_stop = nsurf, trace ray to this surface
% raypos_start = [x y z] to launch ray (m)
% raycos_start = [l m n] to launch ray
% 
% outputs:
% vigcode
% raypos = (x,y,z) (m) (assumes ZEMAX lens units = mm)
% raycos = ray angles (cosx, cosy, cosz)
% surfnorm = surface normal direction (cosx,cosy,coxz) 
% rayI = intensity scalar

global MM;

if nargin == 0,
   error('usage: [vigcode, raypos, raycos, surfnorm, rayI] = z_gettracedirect(zchan, wave, mode, nsurf_start, nsurf_stop, raypos_start, raycos_start)');
end

cmdstr = ['GetTraceDirect,'...
      num2str(wave) ','...
      num2str(mode) ','...
      num2str(nsurf_start) ','...
      num2str(nsurf_stop) ','...
      sprintf('%d,%d,%d,%d,%d,%d',raypos_start/MM,raycos_start)];

rays = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf])';

% rays(1) = error
% rays(2) = vigcode
% rays(3:5) = x,y,z
% rays(6:8) = l,m,n (direction cosines after refraction)
% rays(9:11) = l2,m2,n2 (surface intercept direction normals)
% rays(12) = intensity
status = rays(1);
vigcode = rays(2);

% check status
if status == 0,
    raypos = rays(3:5)*MM; % lens units = MM
    raycos = rays(6:8);
    surfnorm = rays(9:11);
    rayI = rays(12);
else,
    raypos = NaN(3,1);
    raycos = raypos;
    surfnorm = raypos;
    rayI = NaN;

    if status > 0,
        warning(['ray trace error: missed surface ' num2str(status)]);
    else
        warning(['ray trace error: TIR at surface ' num2str(-status)]);
    end
end

return
