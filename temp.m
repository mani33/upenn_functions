function temp(videoFile,outFileName)
% To detect mouse motion.
% Mani Subramaniyan
% 2016-07-01

% Create system objects used for reading video, detecting moving objects,
% and displaying the results.

obj.reader = vision.VideoFileReader(videoFile);
p.FileFormat = 'mp4';
v = VideoWriter(outFileName);
open(v)
i = 0;
while i < (60*30)
    i = i+1;
    frame = obj.reader.step();
   writeVideo(v,frame)  
    if mod(i,30)==0
    disp(['seconds done: ' num2str(round(i/30))]);
    end
end

close(v)



