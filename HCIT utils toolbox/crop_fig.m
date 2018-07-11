
function [] = crop_fig(dname)

  display(['Cropping ',dname '...'])
  F = imread(dname,'png');

  mask = zeros(size(F,1), size(F,2));   % initialize mask
  for ii=1:3,   % for each image plane
    tmp = F(:,:,ii);
    m = tmp~=255;   % create mask
    mask = mask+m;  % add to total mask
  end

  [yy xx] = find(~~mask);
  xmin = min(xx);
  xmax = max(xx);
  ymin = min(yy);
  ymax = max(yy);

  G = F(ymin:ymax,xmin:xmax,:);
  imwrite(G, [dname], 'png')

return

