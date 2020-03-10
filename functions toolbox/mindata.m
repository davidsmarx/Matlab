function varargout = mindata(x,varargin)
% varargout = mindata(x,varargin)
%
% [ap, bp, cp, ...] = filterdata(iuse, a, b, c, ...)
% where [~, iuse] = min(x);
%
% ap = filterdata(iuse, a{:})
%

[~, imin] = min(x);
varargout = cell(size(varargin));
[varargout{:}] = filterdata(imin, varargin{:});
