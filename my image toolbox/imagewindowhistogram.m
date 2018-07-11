function [image_data, h_fig] = imagewindowhistogram(image_input,optionslist)
% [image_data, h_fig] = imagewindowhistogram(image_input,optionslist)
% 
% image_input is a grayscale image array with pixel values 0 to 255
% or it is a filename string for such an image
%
% optionslist is a cell array with any or none of the following strings:
%    'histogram' (default)
%    'profile_x'
%    'profile_y'
%    'surf'
%
% After the image is displayed in a figure, click and drag the mouse to
% draw a rectangle in the image. A histogram of the gray levels for the
% pixels within the rectangle will then be plotted.
%
% outputs:
%   image_data is a matrix
%   h_fig is the handle to the image figure

if isstr(image_input),
    titlestr = pwd2titlestr(image_input);
    image_data = imread(image_input);
else,
    titlestr = inputname(1);
    image_data = image_input;
end

% check for optionslist
if ~exist('optionslist','var'),
    optionslist = {'histogram'}; % default
end

% create the figure for the image and set the callback function for
% mouse down, pass the image to the callback function
h_fig = figure('WindowButtonDownFcn',{@wbDownCallback,image_data,optionslist});

image(image_data), axis image, title(titlestr)

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wbDownCallback(obj,eventdata,imagedata,optionslist)

% get the current point where the mouse down occurred
cpoint = get(gca,'CurrentPoint');
ptstart = cpoint(1,1:2);

% now set the callback for moving the mouse  and pass the image and the 
% current point to the callback
set(gcf,'WindowButtonMotionFcn',{@wbMotionCallback,ptstart,imagedata,optionslist});

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wbMotionCallback(obj,eventdata,ptstart,imagedata,optionslist)

% set up handles for the lines that draw the rectangle
persistent h_line_1 h_line_2 h_line_3 h_line_4;
if isempty(h_line_1),
    h_line_1 = line('erasemode','xor');
    h_line_2 = line('erasemode','xor');
    h_line_3 = line('erasemode','xor');
    h_line_4 = line('erasemode','xor');
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

% draw the rectangle
set(h_line_1,'xdata',[x1 x2],'ydata',[y1 y1]);
set(h_line_2,'xdata',[x1 x1],'ydata',[y1 y2]);
set(h_line_3,'xdata',[x2 x1],'ydata',[y2 y2]);
set(h_line_4,'xdata',[x2 x2],'ydata',[y2 y1]);

% set the callback for when the mouse button is up, and pass the imagedata
% within the rectangle
set(gcf,'WindowButtonUpFcn',{@wbUpCallback,imagedata(y1:y2,x1:x2),optionslist});

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wbUpCallback(obj,eventdata,imageselect,optionslist)
% keep the handle for the histogram figure
persistent h_struct titlestr;

% now that the mouse button is up, end use of the Motion Callback
set(gcf,'WindowButtonMotionFcn',[]);

% set title string
if isempty(titlestr),
    % get title from image
    h_title = get(gca,'title');
    titlestr = get(h_title,'string');
end

if strmatch('histogram',optionslist),
    % produce histogram
    [imhist, graybin] = hist(double(imageselect(:)),256);
    disp(['total pixels: ' num2str(sum(imhist(:)))]);
    if ~isfield(h_struct,'h_hist'), h_struct.h_hist = figure; end
    figure(h_struct.h_hist), bar(graybin,imhist), grid, title(titlestr)
end

if strmatch('profile_x',optionslist),
    % plot intensity profile through maximum point
    [imax, nr, nc] = max2d(double(imageselect));
    if ~isfield(h_struct,'h_profx'), h_struct.h_profx = figure; end
    figure(h_struct.h_profx), plot(imageselect(nr,:)), grid, title([titlestr 'profile x'])    
end

if strmatch('profile_y',optionslist),
    % plot intensity profile through maximum point
    [imax, nr, nc] = max2d(double(imageselect));
    if ~isfield(h_struct,'h_profy'), h_struct.h_profy = figure; end
    figure(h_struct.h_profy), plot(imageselect(:,nc)), grid, title([titlestr 'profile y'])    
end

if strmatch('surf',optionslist),
    % make surf plot of chosen area
    if ~isfield(h_struct,'h_surf'), h_struct.h_surf = figure; end
    figure(h_struct.h_surf), surf(double(imageselect)), title(titlestr)
end

return