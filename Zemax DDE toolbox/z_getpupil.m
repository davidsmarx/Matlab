function pupil = z_getpupil(zchan)
% pupil = z_getpupil(zchan)
%
% pupil is a struct:
%   type
%   value
%   ENPD (entrance pupil diameter)
%   ENPP (entrance pupil position)
%   EXPD (exit pupil diameter)
%   EXPP (exit pupil position)
%   Apodization (struct)
%       type
%       factor
%
% pupil types:
%    0 => entrance pupil diameter
%    1 => image space F#
%    2 => object space NA
%    3 => float by stop size

global MM; % lens units

cmdstr = 'GetPupil';

a = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);

pupil.type = a(1);
pupil.value = a(2);
pupil.ENPD = a(3)*MM; % (entrance pupil diameter)
pupil.ENPP = a(4)*MM; % (entrance pupil position)
pupil.EXPD = a(5)*MM; % (exit pupil diameter)
pupil.EXPP = a(6)*MM; % (exit pupil position)
pupil.Apodization.type = a(7);
pupil.Apodization.factor = a(8);

% pupil value units depends on the type of pupil
switch pupil.type
    case {0, 3},
        pupil.value = pupil.value*MM;
    otherwise,
end

return


