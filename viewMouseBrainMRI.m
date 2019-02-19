function varargout = viewMouseBrainMRI(varargin)
% VIEWMOUSEBRAINMRI MATLAB code for viewMouseBrainMRI.fig
%      VIEWMOUSEBRAINMRI, by itself, creates a new VIEWMOUSEBRAINMRI or raises the existing
%      singleton*.
%
%      H = VIEWMOUSEBRAINMRI returns the handle to a new VIEWMOUSEBRAINMRI or the handle to
%      the existing singleton*.
%
%      VIEWMOUSEBRAINMRI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEWMOUSEBRAINMRI.M with the given input arguments.
%
%      VIEWMOUSEBRAINMRI('Property','Value',...) creates a new VIEWMOUSEBRAINMRI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before viewMouseBrainMRI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to viewMouseBrainMRI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help viewMouseBrainMRI

% Last Modified by GUIDE v2.5 05-Jun-2015 17:12:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @viewMouseBrainMRI_OpeningFcn, ...
                   'gui_OutputFcn',  @viewMouseBrainMRI_OutputFcn, ...
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


% --- Executes just before viewMouseBrainMRI is made visible.
function viewMouseBrainMRI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to viewMouseBrainMRI (see VARARGIN)

% Choose default command line output for viewMouseBrainMRI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes viewMouseBrainMRI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = viewMouseBrainMRI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load_file_list.
function load_file_list_Callback(hObject, eventdata, handles)
% hObject    handle to load_file_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.imageFiles,handles.imageDir] = uigetfile('*.*','Choose slice files','multiselect','on');
nF = length(handles.imageFiles);
for i = 1:nF
    handles.imagePath{i} = fullfile(handles.imageDir,handles.imageFiles{i});    
end
set(handles.loaded_file_list,'String',handles.imageFiles);
guidata(hObject,handles)



% --- Executes on button press in slicer_line.
function slicer_line_Callback(hObject, eventdata, handles)
% hObject    handle to slicer_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[~,c] = size(handles.slicer_im);
[x,y] = ginput(2);
X = [x ones(size(x))];
B = regress(y,X);
xi = [(1:c)' ones(c,1)];
yi = xi*B; 
xx = xi(:,1);
yy = round(yi);

axes(handles.slicer_axes)
imagesc(handles.slicer_im);
axis image; axis xy;
hold on
plot(xx,yy,'y.')
hold off

% % Check if the indices are within bound
% % xinbound = (min(xx) >= 1) && (max(xx) <= c);
% yinbound = (min(yy) >= 1) && (max(yy) <= r);
% if ~yinbound
%     Y = [y ones(size(y))];
%     B = regress(x,Y);
%     yi = [(1:r)' ones(r,1)];
%     xi = yi*B;
%     xx = round(xi);
%     yy = yi(:,1);
%     pind = find(yy>0);
%     xx = xx(pind);
%     yy = yy(pind);
% end
nF = length(handles.imageFiles);
handles.sliced_im = [];


for i = 1:nF
    im = imread(handles.imagePath{i});
    for j = 1:length(yy);
        handles.sliced_im(i,j) = im(yy(j),xx(j));
    end
end
[r,c] = size(handles.sliced_im);
xi = sqrt(1+B(1)^2)*(0:c-1);
yi = (0:r-1);

axes(handles.sliced_image_axes)

imagesc(xi,yi,handles.sliced_im)
axis image; axis xy; colormap gray

function vox_x_Callback(hObject, eventdata, handles)
% hObject    handle to vox_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vox_x as text
%        str2double(get(hObject,'String')) returns contents of vox_x as a double
handles.vx = str2double(get(hObject,'String'));
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function vox_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vox_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vox_y_Callback(hObject, eventdata, handles)
% hObject    handle to vox_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vox_y as text
%        str2double(get(hObject,'String')) returns contents of vox_y as a double
handles.vy = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function vox_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vox_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vox_z_Callback(hObject, eventdata, handles)
% hObject    handle to vox_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vox_z as text
%        str2double(get(hObject,'String')) returns contents of vox_z as a double
handles.vz = str2double(get(hObject,'String'));
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function vox_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vox_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function distance_measured_Callback(hObject, eventdata, handles)
% hObject    handle to distance_measured (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of distance_measured as text
%        str2double(get(hObject,'String')) returns contents of distance_measured as a double


% --- Executes during object creation, after setting all properties.
function distance_measured_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distance_measured (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in measure_dist.
function measure_dist_Callback(hObject, eventdata, handles)
% hObject    handle to measure_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[x,y] = ginput(2);
xp = abs(diff(x));
yp = abs(diff(y));
xd = xp * handles.vx;
yd = yp * handles.vy;
set(handles.distance_measured,'String',sqrt(xd^2 + yd^2))


% --- Executes on selection change in loaded_file_list.
function loaded_file_list_Callback(hObject, eventdata, handles)
% hObject    handle to loaded_file_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns loaded_file_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from loaded_file_list
 contents = cellstr(get(hObject,'String')) ;
 handles.sliceFile = contents{get(hObject,'Value')};

axes(handles.slicer_axes)
handles.slicer_im = imread(fullfile(handles.imageDir,handles.sliceFile));
imagesc(handles.slicer_im)
axis image
axis xy
colormap gray
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function loaded_file_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loaded_file_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
