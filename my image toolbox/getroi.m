function [roi, x, y, h_line, h_fig] = getroi(I,varargin)
% [roi, x, y, h_line, h_fig] = getroi(I,options)
%
% I can be a filename, image matrix, or axes handle
%
% interactively select a region of interest within the image I
% makes use of the Matlab function getrect
% x = [x1, x2], y = [y1, y2], and roi = I(y1:y2,x1:x2)
% On error, x, y, and roi are returned as empty matrices
% 
% options:
% 'linecolor' 'r','b',etc. sets the color of the box (default = blue)
%
% See also GETRECT IMCROP GETROIFIXEDRECT

% choose the region of interest interactively

if isstr(I), % input is a filename
    I = imread(I); 
    h_fig = figure; imshow(I);
    h_axes = gca;
    
elseif ishandle(I), % input is a handle, check if figure or axes
    switch get(I,'type'),
        case 'figure'
            h_fig = I;
            figure(h_fig);
            h_axes = gca;
        case 'axes'
            h_axes = I;
            h_fig = get(I,'parent');
        case 'image'
            h_axes = get(I,'parent');
            h_fig = get(h_axes,'parent');
        otherwise,
            error(['handle type ' get(I,'type') ' not supported']);
    end
    I = getimage(h_fig);
    
else, % input must be an image
    h_fig = figure; imshow(I)
    h_axes = gca;
end
axes(h_axes); % set current axes

% evaluate input options
structOptions = GetOptions(varargin{:});

rect = getrect(h_axes); % rect = [xmin  ymin width height]
x1 = ceil(rect(1)); x2 = floor(rect(1)+rect(3)-1);
y1 = ceil(rect(2)); y2 = floor(rect(2)+rect(4)-1);

% draw the chosen rectangle in the figure
px = [x1 x1 x2 x2 x1];
py = [y1 y2 y2 y1 y1];
h_line = line(px,py,'color',structOptions.linecolor);

% return the region of interest
% test for intensity or RGB image
% test that x and y are within bounds
Ndim = ndims(I);
switch Ndim
    case 3, % RGB image
        [Ny, Nx, N3] = size(I);
        if y1<1 | y2>Ny | x1<1 | x2>Nx, [roi, x, y] = ReturnError; return, end
        roi = I(y1:y2,x1:x2,:);
    case 2,
        [Ny, Nx] = size(I);
        if y1<1 | y2>Ny | x1<1 | x2>Nx, [roi, x, y] = ReturnError; return, end
        roi = I(y1:y2,x1:x2);
    otherwise,
        warning('dimensions error input image');
        [roi, x, y] = ReturnError;
        return
end

x = [x1, x2];
y = [y1, y2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function structOptions = GetOptions(varargin);

% default values:
structOptions.linecolor = 'b';

for ii = 1:length(varargin),
    switch lower(varargin{ii}),
        case 'linecolor',
            structOptions.linecolor = varargin{ii+1};
        otherwise,
            
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [roi, x, y] = ReturnError

roi = [];
x = [];
y = [];

