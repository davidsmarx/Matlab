function [efl, pwfn, rwfn, pima, pmag] = z_getfirstorderdata(zchan)
% [efl, pwfn, rwfn, pima, pmag] = z_getfirstorderdata(zchan)
%
% efl = effective focal length
% pwfn = paraxial working F#
% rwfn = real working F#
% pima = paraxial image height
% pmag = paraxial magnification

global MM;

cmdstr = 'GetFirst';

a = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);

efl = a(1)*MM;
pwfn = a(2);
rwfn = a(3);
pima = a(4)*MM;
pmag = a(5);

return


