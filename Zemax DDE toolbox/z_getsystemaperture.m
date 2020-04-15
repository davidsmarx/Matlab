function sAp = z_getsystemaperture(zchan)
% sAp = z_getsystemaperture(zchan)
%
% sAp.type
% sAp.stopsurf
% sAp.aperture (the value)

U = CConstants;

try
cmdstr = ['GetSystemAper'];
a = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);

switch a(1)
    case 0,
        sAp.type = 'Entrance Pupil Diameter';
        aunits = U.MM;
    case 1,
        sAp.type = 'Image Space F/#';
        aunits = 1;
    case 2,
        sAp.type = 'Object Space NA';
        aunits = 1;
    case 3,
        sAp.type = 'Float by Stop Size';
        aunits = U.MM;
    case 4,
        sAp.type = 'Paraxial Working F/#';
        aunits = 1;
    case 5,
        sAp.type = 'Object Cone Angle';
        aunits = U.P; % half-angle in degrees
    otherwise,
        error('unknown aperture type');
end

sAp.stopsurf = a(2); % surface #

sAp.aperture = a(3)*aunits; % units depend on type

catch,
   disp('ERROR in z_getconfig!');
   keyboard;
end

return


