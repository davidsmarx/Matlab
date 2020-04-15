function retval = z_setsystemaperture(zchan, sAp)
% retval = z_setsystemaperture(zchan, sAp)
%
% sAp.type
% sAp.stopsurf (surface #)
% sAp.aperture (the value) (units depend on type)
%
% not sure how to interpret the retval

U = CConstants;

try

switch sAp.type
    case 'Entrance Pupil Diameter',
        a(1) = 0;
        aunits = U.MM;
    case 'Image Space F/#',
        a(1) = 1;
        aunits = 1;
    case 'Object Space NA',
        a(1) = 2;
        aunits = 1;
    case 'Float by Stop Size',
        a(1) = 3;
        aunits = U.MM;
    case 'Paraxial Working F/#',
        a(1) = 4;
        aunits = 1;
    case 'Object Cone Angle',
        a(1) = 5;
        aunits = U.P; % half-angle in degrees
    otherwise,
        error('unknown aperture type');
end

a(2) = sAp.stopsurf; % surface #

a(3) = sAp.aperture/aunits;

cmdstr = sprintf('SetSystemAper, %d, %d, %d',a);
retval = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);

catch,
   disp('ERROR in z_getconfig!');
   keyboard;
end

return


