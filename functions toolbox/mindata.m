function varargout = mindata(x,varargin)
% varargout = mindata(x,varargin)
%
% [ap, bp, cp, ...] = mindata(x, a, b, c, ...)
% where [~, imin] = min(x);
% ap = a(imin), bp = b(imin, ...
%

[~, imin] = min(x);
varargout = cell(size(varargin));
[varargout{:}] = filterdata(imin, varargin{:});
