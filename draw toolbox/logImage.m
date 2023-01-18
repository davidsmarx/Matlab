function imgout = logImage(img, varargin)
% log image ds9 style
%
% from http://ds9.si.edu/doc/ref/how.html
%
% CheckOption('alpha', 1000, varargin{:});

alpha = CheckOption('alpha', 1000, varargin{:});
imgout = log(alpha*img+1)./log(alpha);

end

