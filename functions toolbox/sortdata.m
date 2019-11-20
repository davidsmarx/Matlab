function varargout = sortdata(cellA,varargin)
% varargout = sortdata(A,varargin)
%
% [a, b, c, ...] = sortdata({A, B, C, ...}, varargin)
% [a, b, c, ..., isort] = sortdata(...)
%
% sort A:
% [a, isort] = sort(a, ...), and apply sort order
% b = B(isort); c = C(isort); ...
%
% A must be 1-d
% if B, C, etc. ndims >= 2, sort B(isort,:)
%
% varargin are optional arguments passed to sort(a, ...)
%

sortdim = CheckOption('sortdim', 1, varargin{:}); % only know how to do 1 for now

[varargout{1}, isort] = sort(cellA{1}, varargin{:});

nsort = min(nargout, length(cellA));
for ii = 2:nsort,
    if isvector(cellA{ii}),
        varargout{ii} = cellA{ii}(isort);
    else
        varargout{ii} = ndimsort(cellA{ii}, isort, sortdim);
    end
end

if nargout == nsort + 1,
    varargout{nargout} = isort;
end


function Bsort = ndimsort(B, isort, sortdim)
% only know sortdim = 1 for now
sizeB = size(B);
%ndimsB = ndims(B);

Br = reshape(B, [sizeB(1) prod(sizeB(2:end))]);
Brsort = Br(isort, :);
Bsort = reshape(Brsort, sizeB);




