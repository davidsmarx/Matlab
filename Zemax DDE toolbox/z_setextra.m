function valout = z_setextra(zchan,surf,collist,vallist)
% values = z_setextra(zchan,surf,collist,vallist)
%
% surf = surface struct, or surface number

switch class(surf)
    case 'struct',
        surfnum = surf.nsurf;
    case 'double',
        surfnum = surf;
    otherwise,
        error('unknown surf');
end

if length(collist(:)) ~= length(vallist(:)),
    error('number of columns does not equal number of values');
end

valout = zeros(size(vallist));
for ii = 1:length(collist(:)),
    cmdstr = ['SetExtra, ' num2str(surfnum) ', ' num2str(collist(ii)) ', ' num2str(vallist(ii),'%e')];
    valtmp = ddereq( zchan, cmdstr, [1 1]);
    valout(ii) = str2double(valtmp);
end

return

