function [image_data, hnim] = imagewindowmicrometer(image_input,optionsstruct)
% [image_data, h_figure] = imagewindowmicrometer(image_input,optionsstruct)
% 
% image_input is a:
%     black & white image array with pixel values 0 to 255, or
%     a filename string for such an image, or
%     a handle to a figure containing an image
%
% optionsstruct is a struct with any or none of the following fields:
%     dxdy = row vector of pixel size [dx, dy] in the x and y-direction. If
%     not input, the result is in pixel units.
%     units = string describing units of dx and dy
%     titlestr = string to use for the figure title
%
% outputs:
%     hnim = handle to the image figure
%     image_data = matrix of image data
%
% After the image is displayed in a figure, click and drag the mouse to
% measure the  distance between two points in the image

% check for optionslist
if ~exist('optionsstruct','var'),
    optionsstruct = struct; % default: no options
end

% check variable type of input image
if isstr(image_input),
    titlestrtmp = pwd2titlestr(image_input);
    image_data = imread(image_input);
    hnim = createimagefigure(image_data);
elseif ishandle(image_input),
    hnim = image_input;
    titlestrtmp = ' ';
else,
    titlestrtmp = inputname(1);
    image_data = image_input;
    hnim = createimagefigure(image_data);
end
% check for titlestr in input options
if ~isfield(optionsstruct,'titlestr'), optionsstruct.titlestr = titlestrtmp; end

% set the figure properties for the image and set the callback function for
% mouse down, pass the image to the callback function
figure(hnim); title(optionsstruct.titlestr);
set(hnim,'pointer','crosshair');
set(hnim,'WindowButtonDownFcn',{@wbDownCallback,optionsstruct});

return

function h_figure = createimagefigure(image_data,titlestr)

h_figure = figure;
image(image_data), axis image

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wbDownCallback(obj,eventdata,optionsstruct)

% get the current point where the mouse down occurred
cpoint = get(gca,'CurrentPoint');
ptstart = cpoint(1,1:2);

% now set the callback for moving the mouse and pass the image and the 
% current point to the callback
set(gcf,'WindowButtonMotionFcn',{@wbMotionCallback,ptstart,optionsstruct});

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wbMotionCallback(obj,eventdata,ptstart,optionsstruct)

% set up a handle for the line drawn from start to current mouse position
persistent h_line;
if isempty(h_line),
    h_line = line('erasemode','xor');
end

% set starting point from the mouse down
x1 = round(ptstart(1));  y1 = round(ptstart(2));
% get the new current point while the mouse is moving
cpoint = get(gca,'CurrentPoint');
x2 = round(cpoint(1,1)); y2 = round(cpoint(1,2));

% check against image limits
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
x1 = ceil(max([x1 xlim(1)])); x2 = floor(min([x2 xlim(2)]));
y1 = ceil(max([y1 ylim(1)])); y2 = floor(min([y2 ylim(2)]));

% draw the line
set(h_line,'xdata',[x1 x2],'ydata',[y1 y2]);

% set the callback for when the mouse button is up, and pass the imagedata
% within the rectangle
set(gcf,'WindowButtonUpFcn',{@wbUpCallback,ptstart,optionsstruct});

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wbUpCallback(obj,eventdata,ptstart,optionsstruct)
% keep the handle for the histogram figure
persistent h_struct;

% now that the mouse button is up, end use of the Motion Callback
set(gcf,'WindowButtonMotionFcn',[]);

% get the mouse point
cpoint = get(gca,'CurrentPoint');
ptstop = cpoint(1,1:2);

% check scaling
if ~isfield(optionsstruct,'dxdy'),
    optionsstruct.dxdy = [1, 1];
    optionsstruct.units = 'pixels';
end
% calculate distance and scale (if dx and dy are set)
cvect = (ptstop-ptstart).*optionsstruct.dxdy;
dist = sqrt(cvect*cvect');
fprintf(['distance equals %f ' optionsstruct.units '\n'],dist);

return