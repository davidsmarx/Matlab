function [roi, roirect] = getroifixedrect(varargin)
% [roi, roirect] = getroifixedrect(image_input,rectsize)
% [roi, roirect] = getroifixedrect(x,y,image_input,rectsize)
% 
% image_input is a grayscale or rgb image array
% or it is a filename string for such an image
%
% x, y are vectors like in image(x,y,image_input)
%
% rectsize is [length_x length_y] in axes units
%
% After the image is displayed in a figure, click and drag the mouse to
% place a rectangle in the image. The rectangle has a fixed size of
% length_x by length_y
%
% outputs:
%   roi is the imaged data within the rectangle
%   roirect = [x1 x2; y1 y2]
%
% See also GETRECT IMCROP GETROI

[x, y, image_data, rectsize, isrgb] = ValidateInput(varargin{:});

% create the figure for the image and set the callback function for
% mouse down, pass the image to the callback function
h_fig = figure;
h_image = imshow(image_data,'XData',x,'YData',y);
h_axes = get(h_image,'parent');
set(h_image,'tag','thisimage');
set(h_axes,'tag','axesimage');

handles = guihandles(h_fig);
guidata(h_fig,handles);

set(h_fig,'WindowButtonDownFcn',{@wbDownCallback,handles,rectsize});

uiwait(h_fig);

handles = guidata(h_fig);
xdata = get(handles.linerect,'xdata');
ydata = get(handles.linerect,'ydata');
x1 = xdata(1); x2 = xdata(3);
y1 = ydata(1); y2 = ydata(3);
if isrgb,
    roi = image_data(y>=y1&y<=y2,x>=x1&x<=x2,:);
else,
    roi = image_data(y>=y1&y<=y2,x>=x1&x<=x2);
end
roirect = [x1 x2; y1 y2];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y, image_data, rectsize, isrgb] = ValidateInput(varargin)

switch length(varargin),
    case 2,
        image_input = varargin{1};
        rectsize = varargin{2};
    case 4,
        x = varargin{1};
        y = varargin{2};
        image_input = varargin{3};
        rectsize = varargin{4};
    otherwise,
        error('invalid # of inputs');
end

if length(rectsize) ~= 2, error('rectsize must be a vector of length 2'); end

% check type of input image
if isa(image_input,'char'), % filename
    image_data = imread(image_input);
else,
    image_data = image_input;
end

switch ndims(image_data),
    case 2,
        [ny, nx] = size(image_data);
        isrgb = 0;
    case 3,
        [ny, nx, n3] = size(image_data);
        if n3 ~= 3, error('invalid input image'); end
        isrgb = 1;
    otherwise,
        error('getroifixedrect: invalid input image');
end

if ~exist('x','var'), x = [1:nx]; end
if ~exist('y','var'), y = [1:ny]; end
if length(y) ~= ny, error('y must be vector of length ny'); end
if length(x) ~= nx, error('x must be vector of length nx'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wbDownCallback(hObject,eventdata,handles,rectsize)

% get the current point where the mouse down occurred
cpoint = get(handles.axesimage,'CurrentPoint');
[x1, x2, y1, y2] = GetPoints(cpoint,rectsize);

% draw the rectangle and add the line handles
h_linerect = line([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'erasemode','xor');
handles.linerect = h_linerect;
guidata(hObject,handles);

% now set the callback for moving the mouse  and pass the image and the 
% current point to the callback
set(hObject,'WindowButtonMotionFcn',{@wbMotionCallback,handles,rectsize});

% set the callback for when the mouse button is up, and pass the imagedata
% within the rectangle
set(hObject,'WindowButtonUpFcn',{@wbUpCallback,handles});

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wbMotionCallback(hObject,eventdata,handles,rectsize)

% get the new current point while the mouse is moving
cpoint = get(handles.axesimage,'CurrentPoint');
[x1, x2, y1, y2] = GetPoints(cpoint,rectsize);

% % check against image limits
% xlim = get(gca,'xlim');
% ylim = get(gca,'ylim');
% x1 = ceil(max([x1 xlim(1)])); x2 = floor(min([x2 xlim(2)]));
% y1 = ceil(max([y1 ylim(1)])); y2 = floor(min([y2 ylim(2)]));

% draw the rectangle
set(handles.linerect,'xdata',[x1 x2 x2 x1 x1],'ydata',[y1 y1 y2 y2 y1]);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x1, x2, y1, y2] = GetPoints(cpoint,rectsize)

ptstart = cpoint(1,1:2);
xc = ptstart(1); yc = ptstart(2);
Nx = rectsize(1); Ny = rectsize(2);
%[-floor(N/2) ceil(N/2)-1]
x1 = round(xc)-floor(Nx/2); x2 = round(xc) + ceil(Nx/2)-1;
y1 = round(yc)-floor(Ny/2); y2 = round(yc) + ceil(Ny/2)-1;

% x1 = round(ptstart(1)); y1 = round(ptstart(2));
% x2 = x1 + rectsize(1) - 1;  y2 = y1 + rectsize(2) - 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wbUpCallback(hObject,eventdata,handles)

% now that the mouse button is up, end use of the Motion Callback
set(hObject,'WindowButtonMotionFcn',[]);

uiresume(hObject);