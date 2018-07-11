function valout = z_getextra(zchan,surf,collist)
% vals = z_getextra(zchan,surf,collist)
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

valout = zeros(size(collist));
for ii = 1:length(collist(:)),
    cmdstr = ['GetExtra, ' num2str(surfnum) ', ' num2str(collist(ii))];
    valtmp = ddereq( zchan, cmdstr, [1 1]);
    valout(ii) = str2double(valtmp);
end

return

