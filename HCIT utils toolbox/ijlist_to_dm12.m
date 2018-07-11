function [idm1, jdm1, idm2, jdm2] = ijlist_to_dm12(ijlist, bDebug)
% [idm1, jdm1, idm2, jdm2] = ijlist_to_dm12(ijlist)
%
% ijlist can be array 1 x 3296
% ijlist can be fits file
% ijlist says which actuators to use
%
% output i,j : dm(i,j) is the actuator

if ischar(ijlist), 
    ijlisttmp = fitsread(PathTranslator(ijlist));
    ijlist = ijlisttmp;
end

% python: i,j = ij % self.nact[idm], ij//self.nact[idm]  % 0-offset
i = mod(ijlist,48); j = floor(ijlist/48);
[idm1, jdm1] = filterdata(j<48, i+1, j+1);
[idm2, jdm2] = filterdata(j>=48, i+1, j-48+1);

if exist('bDebug','var') && ~isempty(bDebug),
    figure,
    plot(idm1, jdm1, 'o'), grid, title('DM1 i j list'), xlabel('i'), ylabel('j')
end

