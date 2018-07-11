function varargout = template_gui(varargin)
% TEMPLATE_GUI M-file for template_gui.fig
%      TEMPLATE_GUI, by itself, creates a new TEMPLATE_GUI or raises the existing
%      singleton*.
%
%      H = TEMPLATE_GUI returns the handle to a new TEMPLATE_GUI or the handle to
%      the existing singleton*.
%
%      TEMPLATE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEMPLATE_GUI.M with the given input arguments.
%
%      TEMPLATE_GUI('Property','Value',...) creates a new TEMPLATE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before template_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to template_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help template_gui

% Last Modified by GUIDE v2.5 25-Aug-2003 11:40:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @template_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @template_gui_OutputFcn, ...
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


% --- Executes just before template_gui is made visible.
function template_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to template_gui (see VARARGIN)

% Choose default command line output for template_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes template_gui wait for user response (see UIRESUME)
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
function varargout = template_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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


% --- Executes on button press in pushbutGo.
function pushbutGo_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutGo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

