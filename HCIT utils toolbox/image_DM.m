function h = image_DM(dm, ijlist, titlestr)
% dm is (m)
% assume 48 x 48

U = CConstants;

DM = zeros(48);
DM(ijlist) = dm(ijlist);

rms_ht = rms(dm(ijlist));

himage = imagesc(DM/U.NM); axis image
colorbar, colorbartitle('WFE (nm)')
title(titlestr)

htrms = text(36,3,['rms = ' num2str(rms(dm(ijlist))/U.NM,'%.1f')] ...
    ,'fontsize',14);
htpv  = text(0,3,['p-v = ' num2str(range(dm(ijlist))/U.NM,'%.1f')] ...
    ,'fontsize',14);


if nargout > 0,
    h = [himage, htrms, htpv];
end