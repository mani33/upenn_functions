function plot_histology(key,side,pc,h,tit,varargin)
% Plot the histology atlas image with electrode location mark onto the
% given axis h
%
args.stim = false; 
args = parseVarArgs(args,varargin{:});

mkr = 'O';
if args.stim
    key.chan_num = -1;
    mkr = 'O';
end
% Get cropped image
[ap,im,xy,ox,oy] = fetchn(his.MouseAtlasHippo(sprintf('side = %u',side)) * his.ElecLesionLoc(key),'ap_loc','im','xy_loc','offset_x','offset_y');
if isempty(im)
    fprintf('No histology exists for mouse %u, so deleting the plot panel\n',key.animal_id)
    delete(h)
    return
end
assert(length(xy)==1,'More than one channel provided')
if nargin < 4
    h = gca;
end
figure(round(ap*100+1000))
imshow(im{:})
hold on
xy = xy{:};
nx = xy(1)-ox;
ny = xy(2)-oy;
plot(nx,ny,'color',pc,'markerfacecolor',pc,'marker',mkr,'markersize',4)
title(tit)
