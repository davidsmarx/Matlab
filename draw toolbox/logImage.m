function imgout = logImage(img, varargin)
% imgout = logImage(img, varargin)
% imgout is log(img) ds9 style
%
% from http://ds9.si.edu/doc/ref/how.html
%
% CheckOption('alpha', 1000, varargin{:});

alpha = CheckOption('alpha', 1000, varargin{:});

% img must be scaled 0 to 1
img_sc = (img - min(img(:)))./range(img(:));
imgout = log10(alpha*img_sc+1)./log10(alpha);

end

