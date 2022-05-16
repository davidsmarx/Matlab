function varargout = maxdata(x,varargin)
% varargout = maxdata(x,varargin)
%
% [ap, bp, cp, ...] = maxdata(x, a, b, c, ...)
% where [~, iuse] = max(x);
%

[~, imax] = max(x);
varargout = cell(size(varargin));
[varargout{:}] = filterdata(imax, varargin{:});
