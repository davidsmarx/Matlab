function z_setsystemaperture(zchan, sAp)
% z_setsystemaperture(zchan, sAp)
%
% sAp.type
% sAp.stopsurf (surface #)
% sAp.aperture (the value)

try

switch sAp.type
    case 'Entrance Pupil Diameter',
        a(1) = 0;
    case 'Image Space F/#',
        a(1) = 1;
    case 'Object Space NA',
        a(1) = 2;
    case 'Float by Stop Size',
        a(1) = 3;
    case 'Paraxial Working F/#',
        a(1) = 4;
    case 'Object Cone Angle',
        a(1) = 5;
    otherwise,
        error('unknown aperture type');
end

a(2) = sAp.stopsurf; % surface #

a(3) = sAp.aperture;

cmdstr = sprintf('SetSystemAper, %d, %d, %d',a);
retval = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);

catch,
   disp('ERROR in z_getconfig!');
   keyboard;
end

return


