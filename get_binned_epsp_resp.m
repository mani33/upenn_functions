function [data,mov_idx] = get_binned_epsp_resp(ch_key,args,varargin)
% data = get_binned_epsp_resp(ch_key,bin_width,varargin)
% MS 2017-08-11

args.invert_time = 0;
args.remove_nrem = 0;
args.remove_sleep = 0;
args.keep_just_nrem = 0;
args.keep_just_sleep = 0;
args.keep_just_motion = 0;

args = parseVarArgs(args,varargin{:});
winstr = sprintf('win = %u',args.pre_event_win);
if double(ch_key.session_start_time) ~= 5333866070
    
    data = struct('raw_t',[],'raw_val',[],'bin_t',[],'bin_val',[]);
    tfac = (1e-6)/60; % to convert to min
    
    if ~isempty(fetch(cont.PreEventMotion(ch_key,winstr)))
        ird = fetch(cstim.FepspSlope(ch_key,'smooth_method_num = 0')*cont.PreEventMotion(winstr),'*');
    else
        ird = fetch(cstim.FepspSlope(ch_key,'smooth_method_num = 0')*cont.PreEventMotionInVid(winstr),'*');
    end
    slopes = -[ird.fepsp_slope];
    % Unfiltered slopes
    data.raw_val_unfilt = slopes;
    traw = double([ird.event_ts]);
    t0 = traw(1);
    mov_idx = [ird.dist_var];
    
    ipi = round((median(diff(traw)))*(1e-6)/60);
    
    if args.invert_time
        t0 = traw(end);
    end
    
    %     % Apply motion threshold
    %     if args.pick_th
    %         args.motion_th = quantile(mov_idx,args.motion_th_quantile);
    %     end
    %     slopes(mov_idx < args.motion_th) = nan;
    nsel = mov_idx >= args.motion_idx_bound(1) & mov_idx <= args.motion_idx_bound(2);
    slopes(~nsel) = nan;
    fprintf('mov: %0.2f -%0.2f\n',args.motion_idx_bound(1),args.motion_idx_bound(2))
    
    
    if args.keep_just_motion
        et = [ird.event_ts];
        % Apply motion based selection of slopes
        d = fetch(cont.MovSegManual(ch_key),'*');
        if isempty(d)
            d = fetch(cont.MovSegInVidManual(ch_key),'*');
        end
        nSeg = length(d.seg_begin);
        % Start with a NaN matrix
        tmp.val = nan(size(et));
        if ~isnan(d.seg_begin(1))
            for iSeg = 1:nSeg
                sel = et >= d.seg_begin(iSeg) & et <= d.seg_end(iSeg);
                tmp.val(sel) = slopes(sel);
            end
        end
        slopes = tmp.val;
    end
    
    data.raw_val = slopes;
    data.raw_t = (traw-t0)*tfac; % in minutes
    
    % Get binned responses
    if args.slope_bw == ipi
        rd.y = slopes;
        rd.t = traw;
        rd.se = zeros(size(slopes));
    else
        %         if args.motion_th <=0
        %             bwstr = sprintf('slope_bw = %0.2f',bin_width);
        %             if args.keep_just_motion
        %                 rd = fetch(cstim.SlopeInMovBinned(ch_key,'smooth_method_num = 0',bwstr),'t','y','se');
        %                 if isempty(rd)
        %                     rd = fetch(cstim.SlopeInMovVidBinned(ch_key,'smooth_method_num = 0',bwstr),'t','y','se');
        %                 end
        %             else
        %                 rd = fetch(cstim.SlopeBinned(ch_key,'smooth_method_num = 0',bwstr),'t','y','se');
        %                 data.se = rd.se;
        %             end
        %         else
        [rd.t,rd.y,data.se] = bin_resp(traw,slopes,args);
        % Threshold motion
        %         end
    end
    if isempty(rd)
        warning('no SlopeBinned data found for given key')
        return
    end
    
    data.bin_val = -rd.y;
    tb = double(rd.t);
    data.bin_t = (tb-t0)*tfac; % in minutes
else
    [data, mov_idx] = get_mouse_141_data(ch_key,args);
end


% Remove non-rem
if args.remove_nrem
    d = fetch(cont.NremSegManual(ch_key),'*');
    d.seg_begin = (d.seg_begin-t0)*tfac;
    d.seg_end = (d.seg_end-t0)*tfac;
    nSeg = length(d.seg_begin);
    if ~isnan(d.seg_begin(1))
        for iSeg = 1:nSeg
            data.bin_val(data.bin_t >= d.seg_begin(iSeg) & data.bin_t <= d.seg_end(iSeg)) = NaN;
        end
    end
end

if args.keep_just_nrem
    d = fetch(cont.NremSegManual(ch_key),'*');
    d.seg_begin = (d.seg_begin-t0)*tfac;
    d.seg_end = (d.seg_end-t0)*tfac;
    nSeg = length(d.seg_begin);
    % Start with a NaN matrix
    tmp.val = nan(size(data.bin_val));
    if ~isnan(d.seg_begin(1))
        for iSeg = 1:nSeg
            sel = data.bin_t >= d.seg_begin(iSeg) & data.bin_t <= d.seg_end(iSeg);
            tmp.val(sel) = data.bin_val(sel);
        end
    end
    data.bin_val = tmp.val;
end



% Remove sleep
if args.remove_sleep
    d = fetch(cont.SleepSegManual(ch_key),'*');
    if count(cont.SleepSegManual(ch_key))==0
        d = fetch(cont.SleepSegManualVid(ch_key),'*');
    end
    d.seg_begin = (d.seg_begin-t0)*tfac;
    d.seg_end = (d.seg_end-t0)*tfac;
    nSeg = length(d.seg_begin);
    if ~isnan(d.seg_begin(1))
        for iSeg = 1:nSeg
            data.bin_val(data.bin_t >= d.seg_begin(iSeg) & data.bin_t <= d.seg_end(iSeg)) = NaN;
        end
    end
end

if args.keep_just_sleep
    d = fetch(cont.SleepSegManual(ch_key),'*');
    d.seg_begin = (d.seg_begin-t0)*tfac;
    d.seg_end = (d.seg_end-t0)*tfac;
    nSeg = length(d.seg_begin);
    % Start with a NaN matrix
    tmp.val = nan(size(data.bin_val));
    if ~isnan(d.seg_begin(1))
        for iSeg = 1:nSeg
            sel = data.bin_t >= d.seg_begin(iSeg) & data.bin_t <= d.seg_end(iSeg);
            tmp.val(sel) = data.bin_val(sel);
        end
    end
    data.bin_val = tmp.val;
end

function  [data, mov_idx] = get_mouse_141_data(ch_key,args)
%%
winstr = sprintf('win = %u',args.pre_event_win);
nBad = 20; % 20 pulses were bad
data = struct('raw_t',[],'raw_val',[],'bin_t',[],'bin_val',[]);
tfac = (1e-6)/60; % to convert to min
% Raw
rd = fetch(cstim.FepspSlope(ch_key,'smooth_method_num = 0')*cont.PreEventMotion(ch_key,winstr),'*');
t1 = double([rd.event_ts]);
t1 = t1-t1(1);
ipi = median(diff(t1))*tfac; % interpulse interval in min
[~,ind] = sort(t1);
rd = rd(ind);

mov_idx1 = [rd.dist_var];
slopes = -[rd.fepsp_slope];

% Apply motion threshold
nsel = mov_idx1 >= args.motion_idx_bound(1)  & mov_idx1 <= args.motion_idx_bound(2);
slopes(~nsel) = nan;

% Binned
% Get binned responses
% bwstr = sprintf('slope_bw = %0.2f',bin_width);
% bd = fetch(cstim.SlopeBinned(ch_key,'smooth_method_num = 0',bwstr),'*');
[bd.t,bd.y,bd.se] = bin_resp(t1,slopes,args);

% get rough number of 20 pulse equivalent bins
nBadBin = ceil(nBad*ipi/args.slope_bw);

filler_raw = nan(1,nBad);
filler_bin = nan(1,nBadBin);
data.raw_val = [slopes(1:(end-nBad)) filler_raw];
data.raw_t = t1;
mov_idx1 = [mov_idx1(1:(end-nBad)) filler_raw];
last_t1 = t1(end);

data.se = [bd.se filler_bin];
data.bin_val = -[bd.y(1:(end-nBadBin)) filler_bin];

data.bin_t = double(bd.t)*tfac; % in minutes

%% For mouse 141 recall exp that had data spread over two sessions
key2 = fetch(acq.Ephys(get_key_from_sess_ts('2017-07-17_15-39-25')));
mch_key = key2;
mch_key.chan_num = ch_key.chan_num;
mch_key.slope_win = ch_key.slope_win;
mch_key.smooth_method_num = ch_key.smooth_method_num;

% Get second session for stitching
sd2 = fetch(cstim.FepspSlope(mch_key)*cont.PreEventMotion(mch_key,winstr),'*');
mov_idx2 = [sd2.dist_var];
t = double([sd2.event_ts]);
t = t-t(1);
t2 = last_t1 + diff(t(1:2)) + t;

slopes2 = -[sd2.fepsp_slope];
% Apply motion filter
nsel = mov_idx2 >= args.motion_idx_bound(1)  & mov_idx2 <= args.motion_idx_bound(2);
slopes2(~nsel) = nan;

data.raw_val = [data.raw_val slopes2];
rt = [data.raw_t t2];
data.raw_t = rt*tfac;

% Binned time
[bd2.t,bd2.y,bd2.se] = bin_resp(t2,slopes2,args);

% bd2 = fetch(cstim.SlopeBinned(mch_key,'smooth_method_num = 0',bwstr),'*');
tb = double(bd2.t);
tb2 = tb*tfac; % in minutes
data.bin_t = [data.bin_t tb2];
% Slopes
data.bin_val = [data.bin_val -bd2.y];
data.se = [data.se bd2.se];

mov_idx = [mov_idx1 mov_idx2];

function [t,y,se] = bin_resp(et,slopeVal,args)
% First we get the stimulation event times and bin data based on
% those times

et = double(et);
% Get inter pulse interval
in = median(diff(et)); % in us
% Get number of pulses in the binning interval
bw = args.slope_bw * 60 *(1e6); % to us
pulsesPerBin = round(bw/in);
if pulsesPerBin==1
    disp('Only a single pulse in each bin. So not binning anything')
    t = et;
    y = slopeVal;
    se = zeros(size(t));
    return
end

tot = length(slopeVal);
nSeg = ceil(tot/pulsesPerBin);
y = nan(1,nSeg);
se = y;
t = y;
for iSeg = 1:nSeg
    idx1 = (iSeg-1)*pulsesPerBin + 1;
    idx2 = min(tot,iSeg*pulsesPerBin);
    sel = idx1:idx2;
    selSlopes = slopeVal(sel);
    selTimes = et(sel);
    t1 = selTimes(1);
    t2 = selTimes(end);
    % Get average of slopes within the bin
    y(iSeg) = nanmean(selSlopes);
    se(iSeg) = nanstd(selSlopes)/sqrt(pulsesPerBin);
    t(iSeg) = selTimes(1)+(t2-t1)/2;
end



