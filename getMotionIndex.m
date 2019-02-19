function data = getMotionIndex(videoFile,varargin)
% data = getMotionIndex(videoFile,varargin)
% Vikas Chelur & Mani Subramaniyan
% 2016-08-22
%
args.th = [];
args.frames_to_skip = [];
args.manual_select = false;
args = parseVarArgs(args,varargin{:});

tic
vobj = vision.VideoFileReader(videoFile);
vinfo = info(vobj);
fps = vinfo.VideoFrameRate;
data = struct;
micro_index = 0;
subplot_index = 0;

px = 0;
py = 0;

if isempty(args.th)
   args.th = 3.43;
end
if isempty(args.frames_to_skip)
    args.frames_to_skip = 5;
end
std_th = 3.43;
dispInt = 5; % display interval in minutes

nSkipImageDisplay = floor(fps*60*dispInt);
plotgrid = [4 6];

gw = getGausswin2d(7);
mouse_size = 100; % mouse size in pixels

satisfied = false;
nStdCalc = 30*300;
data.file_name = videoFile;
data.start_time = datestr(now);
data.fps = fps;

%If user would like to select their input
frame = vobj.step();
frame = double(rgb2gray(frame));
data.video_frame_size = [size(frame,1) size(frame,2)];
cf = gpuArray(frame);
cfv = cf(:)';
pixStd = std(cfv);
zcf = ((cf-mean(cfv))/pixStd);

if  args.manual_select
    while ~satisfied
        figure(12221)
        subplot(1,2,1)
        hist(zcf(:),100);
        [std_th,~] = ginput(1);
        szcf = (zcf < std_th);
        szcf = imfilter(szcf,gw);
        [blob_r,blob_c] = find(szcf);
        subplot(1,2,2)
        imshow(szcf)
        hold on
        plot(blob_c,blob_r,'y.')
        redoit = input('Is the threshold ok? If yes, press ENTER, if no, press any other key','s');
        if isempty(redoit)
            satisfied = true;
        end
    end
end

iFrame = 1; % Not zero because we already read one frame

while ~isDone(vobj)
    iFrame = iFrame + 1;
    if mod(iFrame,args.frames_to_skip)==0
        micro_index = micro_index + 1;
        cf = gpuArray(double(rgb2gray(frame)));
        cfv = cf(:)';
        if  mod(iFrame,nStdCalc)==0
            pixStd = std(cfv);
        end
        % z-score the frame
        zcf = (cf-mean(cfv))/pixStd;
        
        %Calculate centerpoint of mouse
        szcf = zcf < std_th;
        szcf = gather(imfilter(szcf,gw));
        [blob_r,blob_c] = find(szcf);
        blob_dist = sqrt((blob_c - px).^2 + (blob_r - py).^2);
        sel = blob_dist < mouse_size;
        blob_c = blob_c(sel);
        blob_r = blob_r(sel);
        cx = median(blob_c);
        cy = median(blob_r);
        
        %Show a visual representation of the mouse
        if mod(iFrame, nSkipImageDisplay)==0 && subplot_index <= prod(plotgrid)
            figure(2)
            subplot_index = subplot_index + 1;
            subplot(plotgrid(1), plotgrid(2), subplot_index)
            imshow(frame)
            hold on;
            plot(blob_c,blob_r,'y.')
            plot(cx,cy,'r*')
        end
        
        %Store data
        data.dist(micro_index) = sqrt((cx-px)^2 + (cy-py)^2);
        data.median_loc(:,micro_index) = [cx; cy];
        
        px = cx;
        py = cy;
        
        if mod(iFrame,30*60)==0
            disp(iFrame/fps/60)
            toc
        end
    end
    frame = vobj.step();
end
tend = toc;
data.end_time = datestr(now);
release(vobj)
data.time_min = (1:args.frames_to_skip:iFrame)*(1/fps)/60;

% Remove the first point as it is not well defined
data.dist(1) = [];
data.median_loc(:,1) = [];
data.time_min(1) = [];
fprintf('Hours taken: %0.2f\n',tend/60/60)


