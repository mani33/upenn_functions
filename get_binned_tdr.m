function data = get_binned_tdr(ch_key,bin_width,args,varargin)
% data = get_binned_tdr(ch_key,bin_width,varargin)
% MS 2017-08-11

args.invert_time = 0;
args = parseVarArgs(args,varargin{:});

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
    
    % Now get the motion data
    bwstr = sprintf('slope_bw = %0.2f',bin_width);
    
    % Full unbinned data
    fd = fetch(cont.TDratio(ch_key),'*');
    if isempty(fd)
        warning('no TDratio data found for given key from mouse: %u',ch_key.animal_id)
        return
    end
    t_raw = fd.t;
    data.raw_t = (t_raw-t0)*tfac; % in minutes
    data.raw_val = fd.td_ratio;
    
    % Get binned responses
    rd = fetch(cont.TdrBinned(ch_key,bwstr),'*');
    if isempty(rd)
        warning('no TdrBinned data found for given key')
        return
    end
    
    data.se = rd.se;
    data.bin_val = rd.y;
    tb = double(rd.t);
    data.bin_t = (tb-t0)*tfac; % in minutes
    
    % Theta and delta separately
    data.th = rd.th;
    data.se_th = rd.se_th;
    data.del = rd.del;
    data.se_del = rd.se_del;
else
    data = get_mouse_141_data(ch_key,args);
end


function  data = get_mouse_141_data(ch_key,args)
%%
nBad = 20; % 20 pulses were bad
data = struct('raw_t',[],'raw_val',[],'bin_t',[],'bin_val',[]);
tfac = (1e-6)/60; % to convert to min

% Always use the electrical pulse onset time as t0 for comparison
rd = fetch(cstim.FepspSlope(ch_key,'smooth_method_num = 0'),'event_ts');
traw = double([rd.event_ts]);
ipi = median(diff(traw))*tfac;
t0 = traw(1);
if args.invert_time
    t0 = traw(end);
end

nBadBin = ceil(nBad*ipi/args.slope_bw);
% Now get the motion data
bwstr = sprintf('slope_bw = %0.2f',args.slope_bw);

% Full unbinned data
fd = fetch(cont.TDratio(ch_key),'*');
if isempty(fd)
    warning('no TDratio data found for given key')
    return
end
t_raw = fd.t;
data.raw_t = (t_raw-t0)*tfac; % in minutes
data.raw_val = fd.td_ratio;

% Get binned responses
bdata = fetch(cont.TdrBinned(ch_key,bwstr),'*');
if isempty(bdata)
    warning('no TdrBinned data found for given key')
    return
end

data.se = bdata.se;
data.bin_val = bdata.y;
tb = double(bdata.t);
data.bin_t = (tb-t0)*tfac; % in minutes

% Remove last 20 pulse worth data.
nPulseKeep = 100;
tend = (nPulseKeep-1)*ipi;
sel_raw = data.raw_t > tend;
sel_bin = data.bin_t > tend;

data.raw_t(sel_raw) = nan;
data.raw_val(sel_raw) = nan;
data.bin_t(sel_bin) = nan;
data.bin_val(sel_bin) = nan;
data.se(sel_bin) = nan;

% Theta and delta separately
data.th = bdata.th;
data.se_th = bdata.se_th;
data.del = bdata.del;
data.se_del = bdata.se_del;

data.th(sel_bin) = nan;
data.se_th(sel_bin) = nan;
data.del(sel_bin) = nan;
data.se_del(sel_bin) = nan;


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

md2 = fetch(cont.TDratio(mch_key),'*');
if isempty(md2)
    warning('no TDratio data found for given key')
    return
end

t = md2.t;
t2 = last_t1 + (((t-t(1))*tfac) + (nBad+1)*ipi);

data.raw_val = [data.raw_val md2.td_ratio];
data.raw_t = [data.raw_t t2];

% Binned time
% Get binned responses
bd2 = fetch(cont.TdrBinned(mch_key,bwstr),'*');
if isempty(bd2)
    warning('no TdrBinned data found for given key')
    return
end

tb = double(bd2.t);
tb2 = (tb-tb(1))*tfac; % in minutes
bin_t2 = last_bin_t1 + ((nBadBin+1)*bin_width) + tb2;
data.bin_t = [data.bin_t bin_t2];
% Slopes
data.bin_val = [data.bin_val bd2.y];
data.se = [data.se bd2.se];
data.th = [data.th bd2.th];
data.se_th = [data.se_th bd2.se_th];
data.del = [data.del bd2.del];
data.se_del = [data.se_del bd2.se_del];


