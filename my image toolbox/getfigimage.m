function [I, h_image] = getfigimage(handle)
% [I, h_image] = getfigimage(handle)
%
% handle = figure handle or axes handle
% 
% I = the image data (CData) in the figure.
% h_image = the handle to the image
%    If handle is a figure, the image is from the current axes.
%    If handle is an axes, the image is the first image from that axes.
%    If no image is found, I = []; h_image = [];
%

switch get(handle,'type')
    case 'figure'
        h_axes = get(handle,'CurrentAxes');
    case 'axes'
        h_axes = handle;
    otherwise
        error('handle is not a figure or axes handle');
end

% now get the image from axes children
hc = get(h_axes,'children');
h_image = hc( strmatch('image',get(hc,'type')) );
if isempty(h_image), error('no image found'); end
if length(h_image)>1,
    warning('multiple images found, selecting the first');
    h_image = h_image(1);
end
I = get(h_image,'CData');
