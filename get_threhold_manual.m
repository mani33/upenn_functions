function varargout = get_threhold_manual(varargin)
% GET_THREHOLD_MANUAL MATLAB code for get_threhold_manual.fig
%      GET_THREHOLD_MANUAL, by itself, creates a new GET_THREHOLD_MANUAL or raises the existing
%      singleton*.
%
%      H = GET_THREHOLD_MANUAL returns the handle to a new GET_THREHOLD_MANUAL or the handle to
%      the existing singleton*.
%
%      GET_THREHOLD_MANUAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GET_THREHOLD_MANUAL.M with the given input arguments.
%
%      GET_THREHOLD_MANUAL('Property','Value',...) creates a new GET_THREHOLD_MANUAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before get_threhold_manual_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to get_threhold_manual_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help get_threhold_manual

% Last Modified by GUIDE v2.5 22-Aug-2016 20:26:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @get_threhold_manual_OpeningFcn, ...
    'gui_OutputFcn',  @get_threhold_manual_OutputFcn, ...
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


% --- Executes just before get_threhold_manual is made visible.
function get_threhold_manual_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to get_threhold_manual (see VARARGIN)

% Choose default command line output for get_threhold_manual
handles.output = hObject;
set(handles.accept_th,'Value',false)
set(handles.stop_video,'Value',false)
set(handles.pause,'Value',false)
vfile = varargin{1};
handles.key = varargin{2};
set(handles.video_filename,'String',vfile)
handles.videofile = vfile;
handles.video_obj = vision.VideoFileReader(vfile);
% default threshold
th = -3.4;
% Show histogram and plot one frame
handles.current_frame = handles.video_obj.step();
handles.frame_width = size(handles.current_frame,2);
handles.curr_th = th;
set(handles.th_value,'String',num2str(handles.curr_th))
% Update handles structure
guidata(hObject, handles);
update_threshold(th,handles);
axes(handles.hist_display)
hold on
plot([1 1]*handles.curr_th,ylim,'r')
hold off
axes(handles.video_display)
imagesc(handles.current_frame)
% UIWAIT makes get_threhold_manual wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = get_threhold_manual_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1} = handles.accept_th;



function video_filename_Callback(hObject, eventdata, handles)
% hObject    handle to video_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of video_filename as text
%        str2double(get(hObject,'String')) returns contents of video_filename as a double


% --- Executes during object creation, after setting all properties.
function video_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to video_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function update_threshold(th,handles)
frame = handles.current_frame;
cf = double(rgb2gray(frame));
cfv = cf(:)';
pixStd = std(cfv);
zcf = ((cf-mean(cfv))/pixStd);
axes(handles.hist_display)
cla
hist(zcf(:),100);
axis tight
axes(handles.mouse_display)
imshow(frame)
hold on
gw = getGausswin2d(7);
szcf = (zcf < th);
szcf = imfilter(szcf,gw);
[blob_r,blob_c] = find(szcf);
plot(blob_c,blob_r,'y.')


% --- Executes on slider movement.
function th_slider_Callback(hObject, eventdata, handles)
% hObject    handle to th_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
th = get(hObject,'Value');
set(handles.th_value,'String',sprintf('%0.2f',th))
update_threshold(th,handles)
handles.curr_th = th;
guidata(hObject,handles)
axes(handles.hist_display)
hold on
plot([1 1]*handles.curr_th,ylim,'r')

% --- Executes during object creation, after setting all properties.
function th_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to th_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in accept_th.
function accept_th_Callback(hObject, eventdata, handles)
% hObject    handle to accept_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
key = handles.key;
key.th = handles.curr_th;
insert(beh.MotionDetTh,key)



% --- Executes on button press in play_video.
function play_video_Callback(hObject, eventdata, handles)
% hObject    handle to play_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.stop_video,'Value',false)
obj = handles.video_obj;
update_plots = true;
% border_pix_x = 10;
frame = double(obj.step());
gw = getGausswin2d(15);
gate = true;
% frame_width = handles.frame_width;
iFrame = 1;
mouse_size = 100; % mouse size in pixels
axes(handles.mouse_display)
[pcx, pcy] = ginput(1);
nSkip = 5;

while (~obj.isDone()) && ~ get(handles.stop_video,'Value')
    iFrame = iFrame + 1;
    if mod(iFrame,nSkip)==0
        if get(handles.pause,'Value')
            handles.current_frame = frame;
            if (~get(handles.track_button, 'Value')) && gate
                axes(handles.video_display)
                imagesc(frame)
                gate = false;
            end
            if update_plots
                update_threshold(handles.curr_th,handles)
                update_plots = false;
            end
            guidata(hObject,handles)
            pause(1)
        else
            frame = double(obj.step());
            axes(handles.video_display)
            set(handles.frame_num,'String',num2str(iFrame))
            imagesc(frame)
            if get(handles.track_button, 'Value')
                cf = double(rgb2gray(frame));
                cfv = cf(:)';
                pixStd = std(cfv);
                zcf = ((cf-mean(cfv))/pixStd);
                szcf = (zcf < get(handles.th_slider,'Value'));
                szcf = imfilter(szcf,gw);
                [blob_r,blob_c] = find(szcf);
                % Center blob coordinates to previously detected mouse center
                % and take pixels that are only within 150 pixels
                blob_dist = sqrt((blob_c - pcx).^2 + (blob_r - pcy).^2);
                sel = blob_dist < mouse_size;
                blob_c = blob_c(sel);
                blob_r = blob_r(sel);
                hold on
                plot(blob_c,blob_r,'g.')
                cx = median(blob_c);
                cy = median(blob_r);
                plot(cx,cy,'O','markersize',8,'color',[165 42 42]/255,'markerfacecolor',[165 42 42]/255)
                hold off
                pcx = cx;
                pcy = cy;
            end
            update_plots = true;
            drawnow
            
        end
    end
end







% --- Executes on button press in track_button.
function track_button_Callback(hObject, eventdata, handles)
% hObject    handle to track_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of track_button


% --- Executes on button press in pause.
function pause_Callback(hObject, eventdata, handles)
% hObject    handle to pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pause


% --- Executes on button press in stop_video.
function stop_video_Callback(hObject, eventdata, handles)
% hObject    handle to stop_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stop_video
set(handles.pause,'Value',false)


function th_value_Callback(hObject, eventdata, handles)
% hObject    handle to th_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of th_value as text
%        str2double(get(hObject,'String')) returns contents of th_value as a double
th = str2double(get(handles.th_value, 'String'));
update_threshold(th, handles);
set(handles.th_slider, 'Value', th);

% --- Executes during object creation, after setting all properties.
function th_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to th_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in debug_code.
function debug_code_Callback(hObject, eventdata, handles)
% hObject    handle to debug_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over video_filename.
function video_filename_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to video_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = strrep(handles.videofile, '\VT1.mpg', '');
winopen(str)



function frame_num_Callback(hObject, eventdata, handles)
% hObject    handle to frame_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_num as text
%        str2double(get(hObject,'String')) returns contents of frame_num as a double


% --- Executes during object creation, after setting all properties.
function frame_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
