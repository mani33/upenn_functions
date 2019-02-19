function data = get_binned_motion(ch_key,args,varargin)
% data = get_binned_motion(ch_key,bin_width,varargin)
% MS 2017-08-11
% args.motion_idx_quantiles = [0 0.95];
% args.motion_th_quantile = 0.5;
% args.use_pre_event_motion = true;
% args.invert_time = 0;
% args.slope_bw = 5; % sec binning
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
    
    if args.use_pre_event_motion
        [data.raw_val, t_raw] = fetchn(cont.PreEventMotion(ch_key,winstr),'dist_var','event_ts');
        if isempty(data.raw_val)
            [data.raw_val, t_raw] = fetchn(cont.PreEventMotionInVid(ch_key,winstr),'dist_var','event_ts');
        end
        % Filter by motion thresholding
        nsel = data.raw_val >= args.motion_idx_bound(1) & data.raw_val <= args.motion_idx_bound(2);
        data.raw_val(~nsel) = nan;
        
        t_raw = double(t_raw);
        data.raw_t = (t_raw-t0)*tfac; % in minutes
        % Get binned data
        [data.bin_t,data.bin_val,data.se] = bin_resp(data.raw_t,data.raw_val,args);
    else
        
        % Now get the motion data
        bwstr = sprintf('slope_bw = %0.2f',args.slope_bw);
        
        % Full unbinned data
        fd = get_motion(ch_key);
        if isempty(fd)
            warning('no motion related data found for given key')
            return
        end
        t_raw = fd.t;
        
        data.raw_val = fd.dist_var;
        data.raw_t = (t_raw-t0)*tfac; % in minutes
        % Get binned responses
        bdata = get_motion_binned(ch_key,bwstr);
        if isempty(bdata)
            warning('no MotionBinned data found for given key')
            return
        end
        
        data.se = bdata.se;
        data.bin_val = bdata.y;
        tb = double(bdata.t);
        data.bin_t = (tb-t0)*tfac; % in minutes
    end
else
    data = get_mouse_141_data(ch_key,args);
end



function mbdata = get_motion_binned(key,bwstr)

% First check Motion table
mbdata = fetch(cont.MotionBinned(key,bwstr),'*');
if isempty(mbdata)
    mbdata = fetch(cont.MotionInVidBinned(key,bwstr),'*');
    if isempty(mbdata)
        warning('No motion data found')
    end
end


function mdata = get_motion(key)

% First check Motion table
mdata = fetch(cont.Motion(key),'*');
if isempty(mdata)
    mdata = fetch(cont.MotionInVid(key),'*');
    if isempty(mdata)
        warning('No motion data found')
    end
end


function  data = get_mouse_141_data(ch_key,args)
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
nBadBin = ceil(nBad*ipi/args.slope_bw);
if args.invert_time
    t0 = traw(end);
end
if args.use_pre_event_motion
    [data.raw_val, t_raw] = fetchn(cont.PreEventMotion(ch_key,winstr),'dist_var','event_ts');   
     % Filter by motion thresholding
     nsel = data.raw_val >= args.motion_idx_bound(1) & data.raw_val <= args.motion_idx_bound(2);
    data.raw_val(~nsel) = nan;
        
    t_raw = double(t_raw);
    data.raw_t = (t_raw-t0)*tfac; % in minutes
    % Get binned data
    [data.bin_t,data.bin_val,data.se] = bin_resp(data.raw_t,data.raw_val,args);
else
   
    % Now get the motion data
    bwstr = sprintf('slope_bw = %0.2f',args.slope_bw);
    
    % Full unbinned data
    fd = get_motion(ch_key);
    if isempty(fd)
        warning('no motion related data found for given key')
        return
    end
    t_raw = fd.t;
    data.raw_t = (t_raw-t0)*tfac; % in minutes
    data.raw_val = fd.dist_var;
    
    % Get binned responses
    bdata = get_motion_binned(ch_key,bwstr);
    if isempty(bdata)
        warning('no MotionBinned data found for given key')
        return
    end
    
    data.se = bdata.se;
    data.bin_val = bdata.y;
    tb = double(bdata.t);
    data.bin_t = (tb-t0)*tfac; % in minutes
end

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
if args.use_pre_event_motion
    [md2.dist_var, md2.t] = fetchn(cont.PreEventMotion(mch_key,winstr),'dist_var','event_ts');    
    % Filter by motion thresholding
    nsel = md2.dist_var >= args.motion_idx_bound(1) & md2.dist_var <= args.motion_idx_bound(2);
    md2.dist_var(~nsel) = nan;
        
    md2.t = double(md2.t);
    md2.t = (md2.t-t0)*tfac; % in minutes
    % Get binned data
    [bd2.t,bd2.y,bd2.se] = bin_resp(md2.t,md2.dist_var,args);
else
    md2 = get_motion(mch_key);
    if isempty(md2)
        warning('no motion related data found for given key')
        return
    end
    
    % Binned time
    % Get binned responses
    bd2 = get_motion_binned(mch_key,bwstr);
    if isempty(bd2)
        warning('no MotionBinned data found for given key')
        return
    end
end
% Raw data
t = md2.t;
t2 = last_t1 + diff(t(1:2)) + (t-t(1));
data.raw_val = [data.raw_val(:)' md2.dist_var(:)'];
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
in = median(diff(et)); % in us
% Get number of pulses in the binning interval
bw = args.slope_bw;
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



