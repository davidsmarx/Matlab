function varargout = maxdata(x,varargin)
% varargout = maxdata(x,varargin)
%
% [ap, bp, cp, ...] = maxdata(iuse, a, b, c, ...)
% where [~, iuse] = max(x);
%
% ap = maxdata(iuse, a{:})
%

[~, imax] = max(x);
varargout = cell(size(varargin));
[varargout{:}] = filterdata(imax, varargin{:});
