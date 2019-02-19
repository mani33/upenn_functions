function varargout = measureDistanceOnImage(varargin)
% MEASUREDISTANCEONIMAGE MATLAB code for measureDistanceOnImage.fig
%      MEASUREDISTANCEONIMAGE, by itself, creates a new MEASUREDISTANCEONIMAGE or raises the existing
%      singleton*.
%
%      H = MEASUREDISTANCEONIMAGE returns the handle to a new MEASUREDISTANCEONIMAGE or the handle to
%      the existing singleton*.
%
%      MEASUREDISTANCEONIMAGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEASUREDISTANCEONIMAGE.M with the given input arguments.
%
%      MEASUREDISTANCEONIMAGE('Property','Value',...) creates a new MEASUREDISTANCEONIMAGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before measureDistanceOnImage_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to measureDistanceOnImage_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help measureDistanceOnImage

% Last Modified by GUIDE v2.5 27-Jun-2017 09:43:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @measureDistanceOnImage_OpeningFcn, ...
                   'gui_OutputFcn',  @measureDistanceOnImage_OutputFcn, ...
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


% --- Executes just before measureDistanceOnImage is made visible.
function measureDistanceOnImage_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to measureDistanceOnImage (see VARARGIN)

% Choose default command line output for measureDistanceOnImage
handles.output = hObject;
handles.selected_ind = 1;
handles.point_h = [];
handles.pdot_x = nan;
handles.pdot_y = nan;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes measureDistanceOnImage wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = measureDistanceOnImage_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in select_files.
function select_files_Callback(hObject, eventdata, handles)
% hObject    handle to select_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn,pn] = uigetfile('C:\Users\Mani\Desktop\MouseAtlasImages\AtlasHippocampus\*.*','Select Files','','MultiSelect','on');
if ischar(fn)
    fn = {fn};
end
set(handles.file_list,'String',fn);
handles.path = pn;
handles.files = fn;
guidata(hObject,handles)
axes(handles.current_image)
handles.current_im = imread(fullfile(pn,fn{1}));
imagesc(handles.current_im)
axis image off
hold on
plot(handles.pdot_x,handles.pdot_y,'mO','markersize',4,'markerfacecolor','m')

% --- Executes on selection change in file_list.
function file_list_Callback(hObject, eventdata, handles)
% hObject    handle to file_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns file_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from file_list

contents = cellstr(get(hObject,'String'));
sv = get(hObject,'Value');
sf = contents{sv};
handles.selected_file = fullfile(handles.path,sf);
handles.selected_ind = sv;
axes(handles.current_image)
handles.current_im = imread(handles.selected_file);
imagesc(handles.current_im)
axis image off


guidata(hObject,handles);
% --- Executes during object creation, after setting all properties.
function file_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in one_mm.
function one_mm_Callback(hObject, eventdata, handles)
% hObject    handle to one_mm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[x,y] = ginput(2);
dx = abs(diff(x));
dy = abs(diff(y));
d = max([dx dy]);
handles.onemmpixels = d;
guidata(hObject,handles);


% --- Executes on button press in measure_distance.
function measure_distance_Callback(hObject, eventdata, handles)
% hObject    handle to measure_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[x,y] = ginput(2);
s = handles.onemmpixels;
dx = abs(diff(x)/s);
dy = abs(diff(y)/s);
dd = sqrt(dx^2 + dy^2);

set(handles.dx_val,'String',dx)
set(handles.dy_val,'String',dy)
set(handles.abs_distance,'String',dd)

function dx_val_Callback(hObject, eventdata, handles)
% hObject    handle to dx_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx_val as text
%        str2double(get(hObject,'String')) returns contents of dx_val as a double


% --- Executes during object creation, after setting all properties.
function dx_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dy_val_Callback(hObject, eventdata, handles)
% hObject    handle to dy_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dy_val as text
%        str2double(get(hObject,'String')) returns contents of dy_val as a double


% --- Executes during object creation, after setting all properties.
function dy_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dy_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function abs_distance_Callback(hObject, eventdata, handles)
% hObject    handle to abs_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of abs_distance as text
%        str2double(get(hObject,'String')) returns contents of abs_distance as a double


% --- Executes during object creation, after setting all properties.
function abs_distance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to abs_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in previous.
function previous_Callback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = get(handles.file_list,'String');
sf = max(1,handles.selected_ind - 1);
handles.selected_file = fullfile(handles.path,contents{sf});
handles.selected_ind = sf;
axes(handles.current_image)
handles.current_im = imread(handles.selected_file);
hold on
imagesc(handles.current_im)
axis image
axis off
plot(handles.pdot_x,handles.pdot_y,'mO','markersize',4,'markerfacecolor','m')
guidata(hObject,handles);

% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = get(handles.file_list,'String');
sf = min(handles.selected_ind + 1,length(contents));
handles.selected_ind = sf;
handles.selected_file = fullfile(handles.path,contents{sf});
axes(handles.current_image)
handles.current_im = imread(handles.selected_file);
hold on
imagesc(handles.current_im)
axis image
axis off
plot(handles.pdot_x,handles.pdot_y,'mO','markersize',4,'markerfacecolor','m')
guidata(hObject,handles);


% --- Executes on button press in goto_start.
function goto_start_Callback(hObject, eventdata, handles)
% hObject    handle to goto_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[x,y] = ginput(1);
handles.go_x = x;
handles.go_y = y;
guidata(hObject,handles)


function goto_x_Callback(hObject, eventdata, handles)
% hObject    handle to goto_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of goto_x as text
%        str2double(get(hObject,'String')) returns contents of goto_x as a double


% --- Executes during object creation, after setting all properties.
function goto_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to goto_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function goto_y_Callback(hObject, eventdata, handles)
% hObject    handle to goto_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of goto_y as text
%        str2double(get(hObject,'String')) returns contents of goto_y as a double


% --- Executes during object creation, after setting all properties.
function goto_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to goto_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function go_Callback(hObject, eventdata, handles)
% hObject    handle to go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x = str2double(get(handles.goto_x,'String'))*handles.onemmpixels;
y = str2double(get(handles.goto_y,'String'))*handles.onemmpixels;

nx = handles.go_x + x;
ny = handles.go_y + y;
axes(handles.current_image)
hold on
handles.point_h(end+1) = plot(nx,ny,'rO','markersize',4,'markerfacecolor','r');
guidata(hObject,handles)

% --- Executes on button press in clear_points.
function clear_points_Callback(hObject, eventdata, handles)
% hObject    handle to clear_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.point_h)
    delete(handles.point_h)
    handles.point_h = [];
    guidata(hObject,handles)
end
% --- Executes on button press in clear_last_point.
function clear_last_point_Callback(hObject, eventdata, handles)
% hObject    handle to clear_last_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.point_h)
    delete(handles.point_h(end))
    handles.point_h(end) = [];
    guidata(hObject,handles)
end


% --- Executes on button press in put_dot.
function put_dot_Callback(hObject, eventdata, handles)
% hObject    handle to put_dot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.pdot_x,handles.pdot_y] = ginput(1);
axes(handles.current_image)
hold on
handles.point_h(end+1) = plot(handles.pdot_x,handles.pdot_y,'mO','markersize',4,'markerfacecolor','m');
guidata(hObject,handles)


% --- Executes on button press in mark_1_mm.
function mark_1_mm_Callback(hObject, eventdata, handles)
% hObject    handle to mark_1_mm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[~,handles.onemm_y] = ginput(1);
guidata(hObject,handles)

function dv_value_Callback(hObject, eventdata, handles)
% hObject    handle to dv_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dv_value as text
%        str2double(get(hObject,'String')) returns contents of dv_value as a double


% --- Executes during object creation, after setting all properties.
function dv_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dv_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in get_depth_value.
function get_depth_value_Callback(hObject, eventdata, handles)
% hObject    handle to get_depth_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[~,y] = ginput(1);
s = handles.onemmpixels;
vzero = handles.onemm_y - s;
dv = (y-vzero)/s;
set(handles.dv_value,'String',dv)


% --- Executes on button press in clear_image.
function clear_image_Callback(hObject, eventdata, handles)
% hObject    handle to clear_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.current_image)
cla