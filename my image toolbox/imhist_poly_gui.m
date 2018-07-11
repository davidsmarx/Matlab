function varargout = imhist_poly_gui(varargin)
% varargout = imhist_poly_gui(image)
%
% image can be a filename or an image matrix
%
% IMHIST_POLY_GUI M-file for imhist_poly_gui.fig
%      IMHIST_POLY_GUI, by itself, creates a new IMHIST_POLY_GUI or raises the existing
%      singleton*.
%
%      H = IMHIST_POLY_GUI returns the handle to a new IMHIST_POLY_GUI or the handle to
%      the existing singleton*.
%
%      IMHIST_POLY_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMHIST_POLY_GUI.M with the given input arguments.
%
%      IMHIST_POLY_GUI('Property','Value',...) creates a new IMHIST_POLY_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imhist_poly_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imhist_poly_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imhist_poly_gui

% Last Modified by GUIDE v2.5 22-Aug-2003 16:10:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imhist_poly_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @imhist_poly_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before imhist_poly_gui is made visible.
function imhist_poly_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to imhist_poly_gui (see VARARGIN)

% Choose default command line output for imhist_poly_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes imhist_poly_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Load Reference and Test Images
% if no input argument, use file open dialog
if length(varargin) == 0,
    mnuFileOpen_Callback(hObject,eventdata,handles);
else,
    % figure out data type and display image specified in varargin{1}
    if ishandle(varargin{1}),
        Im = getimage(varargin{1});
        titstr = get(varargin{1},'title');
    elseif ischar(varargin{1}),
        % read the image, if it is rgb, convert to grayscale
        Im = imread(varargin{1});
        titstr = varargin{1};
    else % must be an image matrix
        Im = varargin{1};
        titstr = ' ';
    end
    if ndims(Im) == 3, Im = rgb2gray(Im); end
    
    % set current axes, then display image
    axes(handles.axesImage);
    h_image = imshow(Im); title(pwd2titlestr(titstr)) % varargin{1} is 4th arg
    
    % add the image handle to the list
    handles.imageImage = h_image;
    %guidata(gcbo,handles);
    guidata(hObject,handles);
end

% --- Outputs from this function are returned to the command line.
function varargout = imhist_poly_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutSelectROI.
function pushbutSelectROI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutSelectROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'lineRoiPoly'),
    if ishandle(handles.lineRoiPoly), delete(handles.lineRoiPoly); end
end

% use either poly line or rect depending on radio buttons
if get(handles.radbutPolyline,'value'),
    [roi, x, y, h_line, h_fig] = getroi_poly(handles.axesImage,'bwonly');
    % roi contains a logical mask for the polygon region
    % get the image data from the axes handle
    % [I, xdata, ydata] = GetImage(handles.axesImage);
    I = get(handles.imageImage,'CData');
    [cnts, vals] = imhist(I(roi));

else,
    [roi, x, y, h_line, h_fig] = getroi(handles.axesImage);
    % roi contains the image within the roi
    [cnts, vals] = imhist(roi);
end

% add the handle for the line to the handles object
handles.lineRoiPoly = h_line;
guidata(gcbo,handles);

% now get histogram for the mask part of the image
figure, plot(vals,cnts), grid, xlabel('pixel value'), ylabel('counts'),...
    title('Intensity Histogram in Selected ROI','fontsize',14),...
    set(gca,'xlim',[0 max(vals)])

% --------------------------------------------------------------------
function mnuFile_Callback(hObject, eventdata, handles)
% hObject    handle to mnuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuFileOpen_Callback(hObject, eventdata, handles)
% hObject    handle to mnuFileOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%% Get the Image %%%%%%%%%%%%%%%%%%%%
[filename, pathname, fi] = uigetfile({'*.tif; *.jpg'},'Select Image');
if isequal(filename, 0), return, end

% read the image, if it is rgb, convert to grayscale
Im = imread([pathname filename]);
if ndims(Im) == 3, Im = rgb2gray(Im); end

% set current axes, then display image
axes(handles.axesImage);
h_image = imshow(Im); title(pwd2titlestr(filename))

% add the image handle to the list
handles.imageImage = h_image;
%guidata(gcbo,handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function mnuFileExit_Callback(hObject, eventdata, handles)
% hObject    handle to mnuFileExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1);

% --- Executes on button press in radbutPolyline.
function radbutPolyline_Callback(hObject, eventdata, handles)
% hObject    handle to radbutPolyline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radbutPolyline

set(handles.radbutRectangle,'value',0);

% --- Executes on button press in radbutRectangle.
function radbutRectangle_Callback(hObject, eventdata, handles)
% hObject    handle to radbutRectangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radbutRectangle

set(handles.radbutPolyline,'value',0);
