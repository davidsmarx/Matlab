function [x, y, X, Y, R, himage] = imagescgrid(Im, varargin)
% [x, y, X, Y, R, himage] = imagescgrid(Im, dx, dy)
%
% use CreateGrid() to get x,y and then imagesc(x,y,Im)

dx = 1; dy = 1;
if ~isempty(varargin),
    [dx, dy] = deal(varargin{:});
end

[x, y, X, Y, R] = CreateGrid(Im, dx, dy);

himage = imagesc(x, y, Im);
