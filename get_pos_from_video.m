function tmp = get_pos_from_video
% MS 2018-10-02
% For mouse Olive, some one forgot to turn on the video tracking on. So I
% need to extract the position from the video directly.

% This is the file: 'Y:\ephys\raw\2018-07-24_10-19-58\VT1.mpg'
% I copied this file to the local directory
% video_file = 'C:\olive_video.mpg';
video_file = 'C:\VT1.mpg';
% key = fetch(acq.Sessions('session_path like "%2018-07-24_10-19-58%"'));
key = fetch(acq.Sessions('session_path like "%2018-07-24_12-24-03%"'));
 vreader = vision.VideoFileReader(video_file); 
 v = info(vreader);
fps = v.VideoFrameRate;
% Get an estimate of the session duration based on ephys
% recording
[be,en] = fetchn(acq.Ephys(key),'ephys_start_time','ephys_stop_time');
dur = (double(en)-double(be))*(1e-6); % sec
nFrames = round(dur*fps);
eof = false;
rnFrames = 0;
try
while ~eof
   rnFrames = rnFrames + 1;
   [frame,eof] = vreader();
   frame = rgb2gray(frame);
%    imagesc(frame)
   [y,x] = find(frame > 0.5);
   tmp.cx(rnFrames) = median(x);
   tmp.cy(rnFrames) = median(y);
%    hold on
%    plot(tmp.cx,tmp.cy,'r*')
%    drawnow
%    hold off
%    clf
    displayProgress(rnFrames,nFrames)
end
% tget time
tmp.t = linspace(double(be),double(en),rnFrames);
catch
    save('tmpdat','tmp')
end
    
% targFileName = '';
% Mat2NlxVT(targFileName,0,1,[],[1 1 1 0 0 0 0],tmp.t,tmp.cx,tmp.cy)

           
            
