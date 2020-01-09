function [ImCube, hfig, hax] = FitsPath2ImCube(pn, varargin)
% ImCube = FitsPath2ImCube(pn)
%
% simple routine to collect all the fits files in a subdir
% for example, making an image cube out of all the fits files in a
% phase retrieval subdir

% options

%
listfn = dir(PathTranslator([pn '/*.fits']));
Nf = length(listfn);

% get image size
imtmp = fitsread(PathTranslator([pn '/' listfn(1).name]),'image');
finfo = fitsinfo(PathTranslator([pn '/' listfn(1).name]));

ImCube = zeros([Nf size(imtmp)]);
camz = zeros(Nf,1);

ImCube(1,:,:) = imtmp;
camz(1) = FitsGetKeywordVal(finfo.PrimaryData.Keywords,'camz');

for ii = 1:Nf

    imtmp = fitsread(PathTranslator([pn '/' listfn(ii).name]),'image');
    finfo = fitsinfo(PathTranslator([pn '/' listfn(ii).name]));

    ImCube(ii,:,:) = imtmp;
    camz(ii) = FitsGetKeywordVal(finfo.PrimaryData.Keywords,'camz');
end


figure, [hfig, hax, sUserData] = ImageCube(ImCube, camz, 'x', 0, 'y', 0);


