function varargout = sortdata(cellA,varargin)
% varargout = sortdata(A,varargin)
%
% [a, b, c, ...] = sortdata({A, B, C, ...}, varargin)
%
% sort A:
% [a, isort] = sort(a, ...), and apply sort order
% b = B(isort); c = C(isort); ...
%
% varargin are optional arguments passed to sort(a, ...)
%

[varargout{1}, isort] = sort(cellA{1}, varargin{:});

nsort = min(nargout, length(cellA));
for ii = 2:nsort,
    varargout{ii} = cellA{ii}(isort);
end

if nargout == nsort + 1,
    varargout{nargout} = isort;
end

