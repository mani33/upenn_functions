function varargout = unitConversion(varargin)
% UNITCONVERSION MATLAB code for unitConversion.fig
%      UNITCONVERSION, by itself, creates a new UNITCONVERSION or raises the existing
%      singleton*.
%
%      H = UNITCONVERSION returns the handle to a new UNITCONVERSION or the handle to
%      the existing singleton*.
%
%      UNITCONVERSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNITCONVERSION.M with the given input arguments.
%
%      UNITCONVERSION('Property','Value',...) creates a new UNITCONVERSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before unitConversion_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to unitConversion_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help unitConversion

% Last Modified by GUIDE v2.5 24-Sep-2014 08:51:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @unitConversion_OpeningFcn, ...
                   'gui_OutputFcn',  @unitConversion_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before unitConversion is made visible.
function unitConversion_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to unitConversion (see VARARGIN)

% Choose default command line output for unitConversion
handles.output = hObject;
handles.inch_val_prev = 0;
handles.milli_val_prev = 0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes unitConversion wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = unitConversion_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function inches_Callback(hObject, eventdata, handles)
% hObject    handle to inches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inches as text
%        str2double(get(hObject,'String')) returns contents of inches as a double
cv = str2double(get(hObject,'String'));
if cv~=handles.inch_val_prev
    set(handles.millimeters,'String',sprintf('%0.5f',inch2milli(cv)))
end
handles.inch_val_prev = cv;

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function inches_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function millimeters_Callback(hObject, eventdata, handles)
% hObject    handle to millimeters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of millimeters as text
%        str2double(get(hObject,'String')) returns contents of millimeters as a double

cv = str2double(get(hObject,'String'));
if cv~=handles.milli_val_prev
    set(handles.inches,'String',sprintf('%0.5f',milli2inch(cv)))
end
handles.milli_val_prev = cv;


guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function millimeters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to millimeters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function v = inch2milli(c)
v = c * 25.4;

function v = milli2inch(c)
v = c/25.4;
