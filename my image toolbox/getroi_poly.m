function [roi, x, y, h_line, h_fig] = getroi_poly(I,varargin)
% [roi, x, y, h_line, h_fig] = getroi_poly(I,options)
%
% interactive image roi selection using a polygon (uses roipoly)
%
% I is an image (grayscale or rgb), figure handle, axes handle, or file
% name.
% options:
% 'linecolor' 'r','b',etc. sets the color of the box (default = blue)
% 'minsupport' returned roi is minimum size to contain the selected area
% 'bwonly' return roi = only the bw logical mask, otherwise roi = the image
%     within the mask area.
%
% roi = image matrix containing the roi
% x,y points that make the boundary

[I, h_axes, h_fig, structOptions] = ParseInputs(I,varargin{:});

% set current axes and call roipoly
axes(h_axes); 
[bw, x, y] = roipoly;

% draw the poly line
h_line = line(x,y,'color',structOptions.linecolor);

% if 'bwonly' option, set roi to the mask, otherwise get the roi from the image
if strmatch('bwonly',varargin),
    roi = bw;
else,
    % use bw as a mask to get the image roi
    roi = eval([class(I) '( zeros(size(I)) )']);
    % roi has the same size and class as I
    roi(bw) = I(bw);
end

% if 'minsupport' option, copy roi to a matrix of minimum size
if strmatch('minsupport',varargin),
    x1 = floor(min(x)); x2 = ceil(max(x));
    y1 = floor(min(y)); y2 = ceil(max(y));
    roimin = roi(y1:y2,x1:x2);
    roi = roimin;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I, h_axes, h_fig, structOptions] = ParseInputs(I,varargin)

if ishandle(I),
    switch get(I,'type'),
        case 'figure'
            h_fig = I;
            h_axes = get(h_fig,'currentaxes');
        case 'axes',
            h_axes = I;
            h_fig = get(h_axes,'parent');
        otherwise, error('handle type not accepted');
    end
    hh = get(h_axes,'children');
    I = get(hh( strmatch('image',get(hh,'type'))),'CData');
else,
    switch class(I),
        case 'char', % input is a filename
            I = imread(I); 
            h_fig = figure; imshow(I);
            h_axes = gca;
        case {'double', 'uint8', 'uint16'},
            % I must be an image matrix
            h_fig = figure;
            imshow(I);
            h_axes = gca;
        otherwise, error('I is unknown class');
    end % switch class(I)
end % if ishandle(I)

% evaluate input options
% default values:
structOptions.linecolor = 'b';

for ii = 1:length(varargin),
    switch lower(varargin{ii}),
        case 'linecolor',
            structOptions.linecolor = varargin{ii+1};
        otherwise,
            
    end
end
