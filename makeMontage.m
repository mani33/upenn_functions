function varargout = makeMontage(varargin)
% MAKEMONTAGE MATLAB code for makeMontage.fig
%      MAKEMONTAGE, by itself, creates a new MAKEMONTAGE or raises the existing
%      singleton*.
%
%      H = MAKEMONTAGE returns the handle to a new MAKEMONTAGE or the handle to
%      the existing singleton*.
%
%      MAKEMONTAGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAKEMONTAGE.M with the given input arguments.
%
%      MAKEMONTAGE('Property','Value',...) creates a new MAKEMONTAGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before makeMontage_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to makeMontage_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help makeMontage

% Last Modified by GUIDE v2.5 30-Oct-2015 17:13:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @makeMontage_OpeningFcn, ...
                   'gui_OutputFcn',  @makeMontage_OutputFcn, ...
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


% --- Executes just before makeMontage is made visible.
function makeMontage_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to makeMontage (see VARARGIN)

% Choose default command line output for makeMontage
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes makeMontage wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = makeMontage_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in select_image_files.
function select_image_files_Callback(hObject, eventdata, handles)
% hObject    handle to select_image_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn,pn] = uigetfile('*.*','Select Files','','MultiSelect','on');
set(handles.file_list,'String',fn);
handles.path = pn;
handles.files = fn;
guidata(hObject,handles)

% --- Executes on selection change in file_list.
function file_list_Callback(hObject, eventdata, handles)
% hObject    handle to file_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns file_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from file_list

contents = cellstr(get(hObject,'String'));
sf = contents{get(hObject,'Value')};
handles.selected_file = fullfile(handles.path,sf);

axes(handles.current_image)
handles.current_im = imread(handles.selected_file);
imagesc(handles.current_im)
axis image
axis off

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


% --- Executes on button press in crop_image.
function crop_image_Callback(hObject, eventdata, handles)
% hObject    handle to crop_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[c,r] = ginput(2);
r = round(r);
c = round(c);
handles.crop_r = r;
handles.crop_c = c;
cim = handles.current_im(r(1):r(2),c(1):c(2),1:3);
axes(handles.cropped_image)
imagesc(cim);
axis image off
guidata(hObject,handles)


% --- Executes on button press in crop_all_images.
function crop_all_images_Callback(hObject, eventdata, handles)
% hObject    handle to crop_all_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n = length(handles.files);
imc = class(handles.current_im);
 r = handles.crop_r;
    c = handles.crop_c;
    z = size(handles.current_im,3);
handles.cropped_imstack = zeros(diff(r)+1,diff(c)+1,z,imc);
for i = 1:n
    ff = fullfile(handles.path,handles.files{i});
    im = imread(ff);
    cim = im(r(1):r(2),c(1):c(2),:);
    handles.cropped_imstack(:,:,:,i) = cim;
end
guidata(hObject,handles)


% --- Executes on button press in make_cropped_image_montage.
function make_cropped_image_montage_Callback(hObject, eventdata, handles)
% hObject    handle to make_cropped_image_montage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure
montage(handles.cropped_imstack)


% --- Executes on button press in make_original_image_stack.
function make_original_image_stack_Callback(hObject, eventdata, handles)
% hObject    handle to make_original_image_stack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure
n = length(handles.files);
ff = cell(1,n);
for i = 1:n
    ff{i} = fullfile(handles.path,handles.files{i});    
end
montage(ff,'size',[3,2])


% --- Executes on button press in flip_horizontal.
function flip_horizontal_Callback(hObject, eventdata, handles)
% hObject    handle to flip_horizontal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n = size(handles.current_im,3);
for i = 1:n
    handles.current_im(:,:,i) = flipud(handles.current_im(:,:,i));
end
guidata(hObject,handles)
axes(handles.current_image)
imagesc(handles.current_im)
axis image off

% --- Executes on button press in flip_vertical.
function flip_vertical_Callback(hObject, eventdata, handles)
% hObject    handle to flip_vertical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n = size(handles.current_im,3);
for i = 1:n
    handles.current_im(:,:,i) = fliplr(handles.current_im(:,:,i));
end
axes(handles.current_image)
imagesc(handles.current_im)
axis image off
guidata(hObject,handles)

% --- Executes on button press in save_current_im.
function save_current_im_Callback(hObject, eventdata, handles)
% hObject    handle to save_current_im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn,pn] = uiputfile(handles.selected_file,'Select a file name to save');
ff = fullfile(pn,fn);
imwrite(handles.current_im,ff);
