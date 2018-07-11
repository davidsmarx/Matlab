function sAp = z_getsystemaperture(zchan)
% sAp = z_getconfig(zchan)
%
% sAp.type
% sAp.stopsurf
% sAp.aperture (the value)

try
cmdstr = ['GetSystemAper'];
a = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);

switch a(1)
    case 0,
        sAp.type = 'Entrance Pupil Diameter';
    case 1,
        sAp.type = 'Image Space F/#';
    case 2,
        sAp.type = 'Object Space NA';
    case 3,
        sAp.type = 'Float by Stop Size';
    case 4,
        sAp.type = 'Paraxial Working F/#';
    case 5,
        sAp.type = 'Object Cone Angle';
    otherwise,
        error('unknown aperture type');
end

sAp.stopsurf = a(2); % surface #

sAp.aperture = a(3);

catch,
   disp('ERROR in z_getconfig!');
   keyboard;
end

return


