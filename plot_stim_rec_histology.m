function plot_stim_rec_histology(chkeys)
% function plot_stim_rec_histology(keys)

% Get the range of AP locations
% For stimulation electrode positions

mkr = 'O';
pcstim = 'r';
pcrec = 'k';

nKeys = length(chkeys);
aps = fetch(his.ElecLesionLoc(acq.Sessions(chkeys),'chan_num = -1'),'*');
all_stim_locs = [aps.ap_loc];

ustim_ap_locs = fliplr(unique(all_stim_locs));
nuStim = length(ustim_ap_locs);

% Rec
apr = fetch(his.ElecLesionLoc(chkeys,'chan_num >= 0'),'*');
all_rec_ap_locs = [apr.ap_loc];
urec_ap_locs = fliplr(unique(all_rec_ap_locs));
nuRec = length(urec_ap_locs);

rows = max([nuRec nuStim]);
gs = [rows 2];
%% Stim
% First plot the empty images
han = struct;
for iStim = 1:nuStim
    han.stim(iStim) = msubplot(iStim,1,gs);
    im = fetch1(his.MouseAtlasHippo(sprintf('ap_loc = %0.3f and side = 0',ustim_ap_locs(iStim))),'im');
    imagesc(im)
    axis image
    colormap gray
    hold on
end

% Now plot the stim locations
for iKey = 1:nKeys
    chkey = chkeys(iKey);
    [ap,xy,ox,oy] = fetchn(his.MouseAtlasHippo('side = 0') * his.ElecLesionLoc(acq.Sessions(chkey),'chan_num = -1'),'ap_loc','xy_loc','offset_x','offset_y');
    % Get the correct plot axis
    idx = ustim_ap_locs==ap;
    axes(han.stim(idx))
    xy = xy{:};
    nx = xy(1)-ox;
    ny = xy(2)-oy;
    plot(nx,ny,'color',pcstim,'markerfacecolor',pcstim,'marker',mkr,'markersize',4)
    hold on
end

%% Rec
% First plot the empty images

for iRec = 1:nuRec
    han.rec(iRec) = msubplot(iRec,2,gs);
    im = fetch1(his.MouseAtlasHippo(sprintf('ap_loc = %0.3f and side = 1',urec_ap_locs(iRec))),'im');
    imagesc(im)
    axis image
    colormap gray
    hold on
end

% Now plot the stim locations
for iKey = 1:nKeys
    chkey = chkeys(iKey);
    [ap,xy,ox,oy] = fetchn(his.MouseAtlasHippo('side = 1') * his.ElecLesionLoc(chkey),'ap_loc','xy_loc','offset_x','offset_y');
    % Get the correct plot axis
    idx = urec_ap_locs==ap;
    axes(han.rec(idx))
    xy = xy{:};
    nx = xy(1)-ox;
    ny = xy(2)-oy;
    plot(nx,ny,'color',pcrec,'markerfacecolor',pcrec,'marker',mkr,'markersize',4)
    hold on
end

