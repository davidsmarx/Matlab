function varargout = maxdata(x,varargin)
% varargout = maxdata(x,varargin)
%
% [ap, bp, cp, ...] = filterdata(iuse, a, b, c, ...)
% where [~, iuse] = max(x);
%
% ap = filterdata(iuse, a{:})
%

[~, imax] = max(x);
varargout = cell(size(varargin));
[varargout{:}] = filterdata(imax, varargin{:});
