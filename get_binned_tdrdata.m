function [data,pre_event_tdr_raw] = get_binned_tdrdata(ch_key,args,varargin)
% data = get_binned_tdrdata(ch_key,bin_width,args,varargin)
% MS 2018-05-06
%
args = parseVarArgs(args,varargin{:});
winstr = sprintf('win = %u',args.pre_event_win);
if double(ch_key.session_start_time) ~= 5333866070
    data = struct('raw_t',[],'raw_val',[],'bin_t',[],'bin_val',[]);
    tfac = (1e-6)/60; % to convert to min
    
    % Always use the electrical pulse onset time as t0 for comparison
    rd = fetch(cstim.FepspSlope(ch_key,'smooth_method_num = 0'),'event_ts');
    traw = double([rd.event_ts]);
    t0 = traw(1);
    if args.invert_time
        t0 = traw(end);
    end
    
    % For plotting purpose, just get pre-event tdr
    pre_event_tdr_raw = fetchn(cont.PreEventTdr(ch_key,winstr),'td_ratio');
    
    if args.use_pre_event_tdr
        [data.raw_val, mov_idx, t_raw] = fetchn(cont.PreEventTdr(ch_key,winstr)*cont.PreEventMotion(ch_key,winstr),'td_ratio','dist_var','event_ts');
        if isempty(data.raw_val)
            [data.raw_val, mov_idx, t_raw] = fetchn(cont.PreEventTdr(ch_key,winstr)*cont.PreEventMotionInVid(ch_key,winstr),'td_ratio','dist_var','event_ts');
        end
        % Filter by motion thresholding
        nsel = mov_idx >= args.motion_idx_bound(1) & mov_idx <= args.motion_idx_bound(2);
        data.raw_val(~nsel) = nan;
        
        t_raw = double(t_raw);
        data.raw_t = (t_raw-t0)*tfac; % in minutes
        % Get binned data
        [data.bin_t,data.bin_val,data.se] = bin_resp(data.raw_t,data.raw_val,args);
    else
        % Full unbinned data
        fd = fetch(cont.TDratio(ch_key),'*');
        if isempty(fd)
            warning('no tdr related data found for given key')
            return
        end
        t_raw = fd.t;
        
        data.raw_val = fd.td_ratio;
        data.raw_t = (t_raw-t0)*tfac; % in minutes
        % Get binned responses
        [data.bin_t,data.bin_val,data.se] = bin_resp(data.raw_t,data.raw_val,args);
               
        tb = double(data.bin_t);
        data.bin_t = (tb-t0)*tfac; % in minutes
    end
else
    [data,pre_event_tdr_raw] = get_mouse_141_data(ch_key,args);
end


function  [data,pre_event_tdr_raw] = get_mouse_141_data(ch_key,args)
%%
nBad = 20; % 20 pulses were bad
data = struct('raw_t',[],'raw_val',[],'bin_t',[],'bin_val',[]);
tfac = (1e-6)/60; % to convert to min
winstr = sprintf('win = %u',args.pre_event_win);
% Always use the electrical pulse onset time as t0 for comparison
rd = fetch(cstim.FepspSlope(ch_key,'smooth_method_num = 0'),'event_ts');
traw = double([rd.event_ts]);
ipi = median(diff(traw))*tfac;
t0 = traw(1);

if args.invert_time
    t0 = traw(end);
end
raw1 = fetchn(cont.PreEventTdr(ch_key,winstr),'td_ratio');
raw1 = raw1(1:100);
if args.use_pre_event_tdr
    [data.raw_val, mov_idx, t_raw] = fetchn(cont.PreEventTdr(ch_key,winstr)*cont.PreEventMotion(ch_key,winstr),'td_ratio','dist_var','event_ts');
    if isempty(data.raw_val)
        [data.raw_val, mov_idx, t_raw] = fetchn(cont.PreEventTdr(ch_key,winstr)*cont.PreEventMotionInVid(ch_key,winstr),'td_ratio','dist_var','event_ts');
    end
    % Filter by motion thresholding
    nsel = mov_idx >= args.motion_idx_bound(1) & mov_idx <= args.motion_idx_bound(2);
    data.raw_val(~nsel) = nan;
    
    t_raw = double(t_raw);
    data.raw_t = (t_raw-t0)*tfac; % in minutes
    % Get binned data
    [data.bin_t,data.bin_val,data.se] = bin_resp(data.raw_t,data.raw_val,args);
else      
    % Full unbinned data
    fd = fetch(cont.TDratio(ch_key),'*');
    if isempty(fd)
        warning('no motion related data found for given key')
        return
    end
    t_raw = fd.t;
    data.raw_t = (t_raw-t0)*tfac; % in minutes
    data.raw_val = fd.td_ratio;
    
    % Get binned responses
    [data.bin_t,data.bin_val,data.se] = bin_resp(data.raw_t,data.raw_val,args);
    tb = double(data.bin_t);
    data.bin_t = (tb-t0)*tfac; % in minutes
end
% Remove last 20 pulse worth data.
nPulseKeep = 100;
tend = (nPulseKeep-1)*ipi;
sel_raw = data.raw_t > tend;
sel_bin = data.bin_t > tend;


% data.raw_t(sel_raw) = nan;
data.raw_val(sel_raw) = nan;
% data.bin_t(sel_bin) = nan;
data.bin_val(sel_bin) = nan;
data.se(sel_bin) = nan;

last_t1 = max(data.raw_t);
last_bin_t1 = max(data.bin_t);
%% For mouse 141 recall exp that had data spread over two sessions
key2 = fetch(acq.Ephys(get_key_from_sess_ts('2017-07-17_15-39-25')));
mch_key = key2;
mch_key.chan_num = ch_key.chan_num;
mch_key.slope_win = ch_key.slope_win;
mch_key.smooth_method_num = ch_key.smooth_method_num;

% Get second session for stitching
% Full unbinned data
raw2 = fetchn(cont.PreEventTdr(mch_key,winstr),'td_ratio');
pre_event_tdr_raw = [raw1(:)' nan(1,20) raw2(:)'];
if args.use_pre_event_motion
    [md2.td_ratio, mov_idx, md2.t] = fetchn(cont.PreEventTdr(mch_key,winstr)*cont.PreEventMotion(mch_key,winstr),'td_ratio','dist_var','event_ts');
    if isempty(md2.td_ratio)
        [md2.td_ratio, mov_idx, md2.t] = fetchn(cont.PreEventTdr(mch_key,winstr)*cont.PreEventMotionInVid(mch_key,winstr),'td_ratio','dist_var','event_ts');
    end
    % Filter by motion thresholding
    nsel = mov_idx >= args.motion_idx_bound(1) & mov_idx <= args.motion_idx_bound(2);
    md2.td_ratio(~nsel) = nan;
        
    md2.t = double(md2.t);
    md2.t = (md2.t-t0)*tfac; % in minutes
    % Get binned data
    [bd2.t,bd2.y,bd2.se] = bin_resp(md2.t,md2.td_ratio,args);
else
    % Full unbinned data
    md2 = fetch(cont.TDratio(mch_key),'*');
    md2.t = double(md2.t);
    md2.t = (md2.t-t0)*tfac; % in minutes
    if isempty(md2)
        warning('no motion related data found for given key')
        return
    end
    
    % Binned time
    % Get binned responses
    
    [bd2.t,bd2.y,bd2.se] = bin_resp(md2.t,md2.td_ratio,args);
    
end
% Raw data
t = md2.t-md2.t(1);
t2 = last_t1 + diff(t(1:2)) + t;
data.raw_val = [data.raw_val(:)' md2.td_ratio(:)'];
data.raw_t = [data.raw_t(:)' t2(:)'];

% Binned data
tb = double(bd2.t);
tb2 = (tb-tb(1)); % in minutes
bin_t2 = last_bin_t1 + diff(tb2(1:2)) + tb2;
data.bin_t = [data.bin_t(:)' bin_t2(:)'];

data.bin_val = [data.bin_val(:)' bd2.y(:)'];
data.se = [data.se(:)' bd2.se(:)'];

function [t,y,se] = bin_resp(et,slopeVal,args)
% First we get the stimulation event times and bin data based on
% those times

et = double(et);
% Get inter pulse interval
in = nanmedian(diff(et)); % in us
% Get number of pulses in the binning interval
bw = args.slope_bw;
pulsesPerBin = round(bw/in);
if pulsesPerBin==1
    disp('Only a single pulse in each bin. So not populating cstim.SlopesBinned table')
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



