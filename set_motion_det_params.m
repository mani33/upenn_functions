function varargout = set_motion_det_params(varargin)
% SET_MOTION_DET_PARAMS MATLAB code for set_motion_det_params.fig
%      SET_MOTION_DET_PARAMS, by itself, creates a new SET_MOTION_DET_PARAMS or raises the existing
%      singleton*.
%
%      H = SET_MOTION_DET_PARAMS returns the handle to a new SET_MOTION_DET_PARAMS or the handle to
%      the existing singleton*.
%
%      SET_MOTION_DET_PARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SET_MOTION_DET_PARAMS.M with the given input arguments.
%
%      SET_MOTION_DET_PARAMS('Property','Value',...) creates a new SET_MOTION_DET_PARAMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before set_motion_det_params_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to set_motion_det_params_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help set_motion_det_params

% Last Modified by GUIDE v2.5 24-Aug-2016 20:28:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @set_motion_det_params_OpeningFcn, ...
    'gui_OutputFcn',  @set_motion_det_params_OutputFcn, ...
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


% --- Executes just before set_motion_det_params is made visible.
function set_motion_det_params_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to set_motion_det_params (see VARARGIN)

% Choose default command line output for set_motion_det_params
handles.output = hObject;
set(handles.accept_params,'Value',false)
set(handles.stop_video,'Value',false)
set(handles.pause,'Value',false)
vfile = varargin{1};
handles.key = varargin{2};
set(handles.video_filename,'String',vfile)
handles.video_file = vfile;
handles.video_obj = vision.VideoFileReader(vfile);
v = info(handles.video_obj);
handles.fps = v.VideoFrameRate;
handles.init_mouse_cx = 0;
handles.init_mouse_cy = 0;
handles = initialize_video_processing(handles);

guidata(hObject,handles)
% UIWAIT makes set_motion_det_params wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = set_motion_det_params_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1} = handles.accept_params;



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


function update_hist(handles)
frame = handles.current_frame;
cf = double(rgb2gray(frame));
cfv = cf(:)';
pixStd = std(cfv);
zcf = ((cf-mean(cfv))/pixStd);
axes(handles.hist_display)
cla
hist(zcf(:),100);
hold on
plot([1 1]*handles.curr_th,ylim,'r')
hold off
axis tight

% --- Executes on slider movement.
function handles = th_slider_Callback(hObject, eventdata, handles)
% hObject    handle to th_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
th = get(hObject,'Value');
set(handles.th_value,'String',sprintf('%0.2f',th))
handles.curr_th = th;
% This is a hack because guidata doesn't seem to save the value
save('curr_th','th')
update_hist(handles)
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function th_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to th_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in accept_params.
function accept_params_Callback(hObject, eventdata, handles)
% hObject    handle to accept_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if logical(get(handles.failed_motion_det,'Value'))
    error('Please create cont.MotionDetIgnore table')
    insert(beh.MotionDetIgnore,handles.key)
else
    load('curr_th')
    handles = guidata(hObject);
    key = handles.key;
    vinfo = info(handles.video_obj);
    key.fps = vinfo.VideoFrameRate;
    key.video_frame_size = vinfo.VideoSize;
    %     key.zscore_th = handles.curr_th;
    key.zscore_th = th;
    key.mouse_radius = str2double(get(handles.mouse_radius,'String'));
    key.smooth_ker_size = get(handles.smooth_ker_size_slider,'Value');
    %     key.frames_to_skip = str2double(get(handles.frames_to_skip,'String'));
    key.video_file = handles.video_file;
    key.strel_size = str2double(get(handles.strel_size,'String'));
    key.init_mcx = handles.init_mouse_cx;
    key.init_mcy = handles.init_mouse_cy;
    insert(cont.MouseVidPosDetParams,key)
end
% --- Executes on button press in play_video.
function play_video_Callback(hObject, eventdata, handles)
% hObject    handle to play_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.stop_video,'Value',false)
obj = handles.video_obj;
update_plots = true;

% frame = double(obj.step());


gate = true;
iFrame = handles.nSkipBegin;
axes(handles.video_display)
px = handles.init_mouse_cx;
py = handles.init_mouse_cy;
st = strel('disk',str2double(get(handles.strel_size,'String')));
% st_rect = strel('rectangle',[20 50]);
data = struct;
j = 0;
while (~obj.isDone()) && ~ get(handles.stop_video,'Value')
    iFrame = iFrame + 1;
    nSkip = str2double(get(handles.frames_to_skip,'String'))+1;
    frame = double(obj.step());
    if  mod(iFrame,nSkip)==0
        gw = getGausswin2d(round(str2double(get(handles.smooth_ker_val,'String'))));
        mouse_radius = str2double(get(handles.mouse_radius,'String'));
        
        if get(handles.pause,'Value')
            st = strel('disk',str2double(get(handles.strel_size,'String')));
            handles.current_frame = frame;
            if gate
                axes(handles.video_display)
                imagesc(frame)
                gate = false;
            end
            if update_plots
                update_hist(handles)
                update_plots = false;
            end
            guidata(hObject,handles)
            pause(1)
        else
            j = j + 1;
            axes(handles.video_display)
            set(handles.frame_num,'String',num2str(iFrame))
            imagesc(frame)
            
            cf = gpuArray((rgb2gray(frame)));
            cfv = cf(:)';
            pixStd = std(cfv);
            zcf = gather(((cf-mean(cfv))/pixStd));
            szcf = zcf < get(handles.th_slider,'Value');%
            % %             szcf = bwareaopen((szcf),1000);
            szcf = (imclose(szcf,st));
            szcf = gather(imfilter(gpuArray(szcf),gw));
            [blob_r,blob_c] = find(szcf);
            data.blob_size(j) = length(blob_r);
            
            if isempty(blob_r)
                % First try simple thresholding
                szcf = gather(zcf < get(handles.th_slider,'Value'));%
                [blob_r,blob_c] = find(szcf);
                mouse_det_size = 0.25*median(data.blob_size);
                mouse_found = length(blob_r)>= mouse_det_size;
                % Adjust threshold until we can detect enough pixels
                if ~mouse_found
                    th_original = get(handles.th_slider,'Value');
                    th_low_lim = min(gather(cfv));
                    th_search_vals = linspace(th_original,th_low_lim,30);
                    
                    for iSearch = 1:length(th_search_vals)
                        curr_th = th_search_vals(iSearch);
                        disp(curr_th)
                        set(handles.temp_th,'String',curr_th)
                        szcf = gather(zcf < curr_th);%
                        [blob_r,blob_c] = find(szcf);
                        if length(blob_r)>= mouse_det_size
                            break
                        end
                    end
                end
            end
            % Center blob coordinates to previously detected mouse center
            % and take pixels that are only within certain pixels
            if isnan(px) && ~isempty(blob_r)
                sel = true(1,length(blob_c));
            else
                blob_dist = sqrt((blob_c - px).^2 + (blob_r - py).^2);
                sel = blob_dist < mouse_radius;
            end
            blob_c = blob_c(sel);
            blob_r = blob_r(sel);
            hold on
            plot(blob_c,blob_r,'g.')
            cx = median(blob_c);
            cy = median(blob_r);
            plot(cx,cy,'O','markersize',8,'color',[165 42 42]/255,'markerfacecolor',[165 42 42]/255)
            data.dist(j) = sqrt((cx-px)^2 + (cy-py)^2);
            plot([px cx],[py cy],'y')
            plot(px,py,'kO','markersize',8)
            
            px = cx;
            py = cy;
            %             disp([px py data.dist(end)])
            hold off
            update_plots = true;
            drawnow
            
            axes(handles.motion_index_display)
            cla
            if length(data.dist)<=750
                plot(data.dist,'r.-')
                ylim([0 100])
            else
                plot(data.dist((end-750):end),'r.-')
                ylim([0 100])
            end
            
        end
    end
end

handles.mouse_cx = px;
handles.mouse_cy = py;
guidata(hObject,handles)

function handles = initialize_video_processing(handles)

% Show histogram and plot one frame


key = handles.key;
be = fetch1(acq.Ephys(key),'ephys_start_time');

ets = double(fetchn(acq.Events(key,'event_ttl in (1,128)'),'event_ts'));
tmin = min(ets);

vinfo = info(handles.video_obj);
fps = vinfo.VideoFrameRate;

% Get motion only from the first electrical pulse to the last
% one

% Get the number of frames to skip to get to the first
% electrical pulse
handles.nSkipBegin = round(((tmin-double(be))*1e-6)*fps);


for iFrame = 1:handles.nSkipBegin
handles.video_obj.step();
end

handles.current_frame = handles.video_obj.step();
cf = double(rgb2gray(handles.current_frame));
cfv = cf(:);
zf = (cfv-mean(cfv))/std(cfv);
handles.curr_th = quantile(zf,0.03);
set(handles.th_slider,'Min',min(zf),'Max',max(zf))
set(handles.th_value,'String',num2str(handles.curr_th))
set(handles.mouse_radius,'String',num2str(100))
set(handles.frames_to_skip,'String',0)
set(handles.smooth_ker_val,'String',13)
set(handles.strel_size,'String',15)
set(handles.failed_motion_det,'Value',false)
% Update handles structure


update_hist(handles);
axes(handles.video_display)
imagesc(handles.current_frame)
axes(handles.motion_index_display)
cla
set(handles.frame_num,'String',num2str(handles.nSkipBegin))


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
release(handles.video_obj) 
for iFrame = 1:handles.nSkipBegin
handles.video_obj.step();
end
set(handles.pause,'Value',false)


function th_value_Callback(hObject, eventdata, handles)
% hObject    handle to th_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of th_value as text
%        str2double(get(hObject,'String')) returns contents of th_value as a double
th = str2double(get(handles.th_value, 'String'));
update_hist(handles);
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


% --- Executes on slider movement.
function smooth_ker_size_slider_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_ker_size_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.smooth_ker_val,'String',round(get(hObject,'Value')))

% --- Executes during object creation, after setting all properties.
function smooth_ker_size_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smooth_ker_size_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function mouse_radius_Callback(hObject, eventdata, handles)
% hObject    handle to mouse_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mouse_radius as text
%        str2double(get(hObject,'String')) returns contents of mouse_radius as a double


% --- Executes during object creation, after setting all properties.
function mouse_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mouse_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frames_to_skip_Callback(hObject, eventdata, handles)
% hObject    handle to frames_to_skip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frames_to_skip as text
%        str2double(get(hObject,'String')) returns contents of frames_to_skip as a double


% --- Executes during object creation, after setting all properties.
function frames_to_skip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frames_to_skip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function smooth_ker_val_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_ker_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smooth_ker_val as text
%        str2double(get(hObject,'String')) returns contents of smooth_ker_val as a double


% --- Executes during object creation, after setting all properties.
function smooth_ker_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smooth_ker_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.init_mouse_cx, handles.init_mouse_cy] = ginput(1);
guidata(hObject,handles)


% --- Executes on button press in clear_dist_plot.
function clear_dist_plot_Callback(hObject, eventdata, handles)
% hObject    handle to clear_dist_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.motion_index_display)
cla



function strel_size_Callback(hObject, eventdata, handles)
% hObject    handle to strel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of strel_size as text
%        str2double(get(hObject,'String')) returns contents of strel_size as a double


% --- Executes during object creation, after setting all properties.
function strel_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in failed_motion_det.
function failed_motion_det_Callback(hObject, eventdata, handles)
% hObject    handle to failed_motion_det (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of failed_motion_det



function temp_th_Callback(hObject, eventdata, handles)
% hObject    handle to temp_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of temp_th as text
%        str2double(get(hObject,'String')) returns contents of temp_th as a double


% --- Executes during object creation, after setting all properties.
function temp_th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to temp_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
