function [xc, yc, ccpeak, ccradius, h_line1, h_line2] = xcorr_roi(I1,I2,varargin)
% [xc, yc, ccpeak, ccradius, h_line1, h_line2] = xcorr_roi(I1,I2,options)
%
% I1, I2 can be filenames, image matrices, or axes handles
% I1 = reference image
% I2 = test image
%
% xc, yc is the translation of I2 from best correlation with I1
% On error, xc = 0, yc = 0
% ccpeak = normalized correlation peak value
% ccradius = radius of correlation peak at 0.5
% h_line1 = handle to line object for roi in I1
% h_line2 = handle to line object for roi in I2
%
% options:
%    'select' = '1', '2', 'both', 'none' is which input is used for
%       interactive selection of a region of interest (default = none)
%    'useedges' option, then edge(I,'sobel') is called for each I, and the
%    resulting edge images are used for cross correlation.
%    options for 'getroi' may also be included as the options list is
%    passed to 'getroi' when it is called.

%check type of inputs I1,I2
I1 = CheckType(I1);
I2 = CheckType(I2);

% parse options
structOptions = GetOptions(varargin{:});

switch structOptions.select,
    case '1',
        % first select roi in reference image (I1)
        [I1roi, x, y, h_line1, h_fig] = getroi(I1,varargin{:});
        if isempty(I1roi), 
            [xc, yc, ccpeak, ccradius, h_line1, h_line2] = ErrorReturn('empty ROI selection');
            return, 
        end
        % get same roi in test image
        I2roi = I2(y(1):y(2),x(1):x(2));
        %h_line2 = DrawLine(h_axes2,x,y,'color',structOptions.linecolor);
    case '2',
        [I2roi, x, y, h_line2, h_fig] = getroi(I2,varargin{:});
        if isempty(I1roi), 
            [xc, yc, ccpeak, ccradius, h_line1, h_line2] = ErrorReturn('empty ROI selection');
            return,
        end
        I1roi = I1(y(1):y(2),x(1):x(2));
        %h_line2 = DrawLine(h_axes2,x,y,'color',structOptions.linecolor);
    case 'both',
        [I1roi, x, y, h_line1, h_fig] = getroi(I1,varargin{:});
        [I2roi, x, y, h_line2, h_fig] = getroi(I2,varargin{:});
        if isempty(I1roi) | isempty(I2roi),
            [xc, yc, ccpeak, ccradius, h_line1, h_line2] = ErrorReturn('empty ROI selection');
            return,
        end
    case 'none',
        % do nothing, pass the input images directly to xcorr
        I1roi = I1;
        I2roi = I2;
    otherwise,
        warning('unrecognized select option');
        xc = 0; yc = 0; return
end

% do the xcorr2, if 'useedges', then call edge routine first
if structOptions.useedges,
    hsobel = fspecial('sobel');
    I1roi = imadd( imfilter(I1roi,hsobel,'replicate'),...
        imfilter(I1roi,hsobel','replicate') );
    I2roi = imadd( imfilter(I2roi,hsobel,'replicate'),...
        imfilter(I2roi,hsobel','replicate') );
end
cc = normxcorr2(I1roi,I2roi);

%define the center--zero offset
[Ncy, Ncx] = size(cc);
x0 = ceil(Ncx/2);
y0 = ceil(Ncy/2);

% find the correlation peak offset from the center
[ccmax, nr, nc] = max2d(cc);
[fm, xm, ym] = findpeak2(cc(nr-1,nc),cc(nr,nc-1),cc(nr,nc),cc(nr,nc+1),cc(nr+1,nc));
xc = nc + xm - x0;
yc = nr + ym - y0;

% determine correlation peak and radius, large peak and small radius mean
% good correlation. If there is more than one peak above the threshold,
% choose the highest peak.
thresh = 0.5;
ccthresh = cc>thresh;
[ccB, ccL] = bwboundaries(ccthresh); % ccb is a cell array of [Nx2] matrices of y,x pairs
if isempty(ccB), 
    [xc, yc, ccpeak, ccradius, h_line1, h_line2] = ErrorReturn('no correlation peak detected'); 
    return, 
end

for ib = 1:length(ccB), % choose boundary region with highest correlation
    ccmax(ib) = max2d(cc(ccL==ib));
end
[ccpeak, ibmax] = max(ccmax);

[ccradius, eqboundary] = FitBoundary(ccB{ibmax},2);

% 
% figure, imshow(ccthresh), axis on, line(eqboundary(:,2),eqboundary(:,1))
% fprintf('radius = %f pixels\n',ccradius);
% fprintf('peak correlation = %f\n',ccpeak);


% figure, image(I1roi), title('ROI #1')
% figure, image(I2roi), title('ROI #2')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eqradius, bout] = FitBoundary(boundary,N)
% boundary = [Mx2] matrix of y,x pairs that describe the boundary
% N = order of fit

% create parametric functions to describe boundary: y = fy(t), x = fx(t)
yy = boundary(:,1); xx = boundary(:,2);
M = length(xx);
t = 2*pi*[0:M-1]'./(M-1);

% represent fy and fx as sums of sines and cosines and calculate MMSE
% coefficients
A = [cos(t*[0:N]) sin(t*[1:N])];
xq = A * (A\xx);
yq = A * (A\yy);
% smoothed boundary is collection of new x and y values
bout = [yq xq];

% estimate perimeter from boundary
delta_sq = diff(bout).^2;    
perimeter = sum(sqrt(sum(delta_sq,2)));
% radius of the circle with the same perimeter
eqradius = perimeter/2/pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I, h_axes] = CheckType(I);

if isstr(I), % filename
    I = imread(I);
    %figure, imshow(I), h_axes = gca;
elseif ishandle(I), % axes handle
    h_axes = I;
    hc = get(h_axes,'children');
    ih = strmatch('image',get(hc,'type'));
    if isempty(ih), error('axes does not contain an image'); end
    I = get(hc(ih),'CData');
else, % must be an image matrix
    %figure, imshow(I), h_axes = gca;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h_line = DrawLine(h_axes,x,y,varargin);

axes(h_axes);
px = [x(1) x(1) x(2) x(2) x(1)];
py = [y(1) y(2) y(2) y(1) y(1)];
h_line = line(px,py,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xc, yc, ccp, ccr, h_line1, h_line2] = ErrorReturn(wstring)

warning(wstring);
xc = 0;
yc = 0;
ccp = 0;
ccr = 0;
h_line1 = [];
h_line2 = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function structOptions = GetOptions(varargin);

% default values:
structOptions.select = 'none';
structOptions.linecolor = 'b';
structOptions.useedges = 0;

for ii = 1:length(varargin),
    switch lower(varargin{ii}),
        case 'select',
            structOptions.select = varargin{ii+1};
        case 'linecolor',
            structOptions.linecolor = varargin{ii+1};    
        case 'useedges',
            structOptions.useedges = 1;
        otherwise,
            
    end
end

