function [nrem_act_ratio, nrem_slopes, act_slopes, nrem_trace, act_trace,t] = get_nrem_act_raio(channel_key,pre_event_win, motion_th_quantile,get_traces)
% Compute nrem over active period slope ratio
% MS 2018-05-07
key = channel_key;
nrem_trace = NaN;
act_trace = NaN;
t = NaN;
nrem_slopes = t;
act_slopes = t;
% Get nonRem segments
winstr = sprintf('win = %u',pre_event_win);
nr = fetch(cont.NremSegManual(key),'*');
sd = fetch((cstim.FepspSlope(key)*cstim.FpRespTrace(key))*cont.PreEventMotion(key,winstr),'*');
if isempty(nr)
    nr = fetch(cont.NremSegInVidManual(key),'*');
    sd = fetch((cstim.FepspSlope(key)*cstim.FpRespTrace(key))*cont.PreEventMotionInVid(key,winstr),'*');
end
if all(isnan(nr.seg_begin))
    nrem_act_ratio = NaN;
    disp('No nrem segment found')
    return
end
sel = ~isnan(nr.seg_begin);
sbegin = nr.seg_begin(sel);
send = nr.seg_end(sel);
nSeg = length(sbegin);

% Get slopes and motion info

% Get all slopes within the nrem segments

et = double([sd.event_ts]);
slopes = -[sd.fepsp_slope];
slopes = slopes(:)';
tmp.nrslopes = cell(1,nSeg);
tmp.nrem_sel = cell(1,nSeg);
for iSeg = 1:nSeg
    s1 = sbegin(iSeg);
    s2 =send(iSeg);
    sel = find(et >= s1 & et < s2);
    tmp.nrslopes{iSeg} = slopes(sel);
    tmp.nrem_sel{iSeg} = sel(:)';
end
nrem_slopes = [tmp.nrslopes{:}];
nrem_mean = mean(nrem_slopes);
nrem_sel = [tmp.nrem_sel{:}];

% Pick slopes for which the pre event motion was above threshold
mov_idx = [sd.dist_var];
mth = quantile(mov_idx,motion_th_quantile);
act_sel = mov_idx > mth;
act_slopes = slopes(act_sel);
if isempty(act_slopes)
    disp('No active region slopes were found')
    nrem_act_ratio = NaN;
else
    act_mean = mean(act_slopes);
    nrem_act_ratio = nrem_mean/act_mean;
end
if get_traces
    % Get waveforms
    nrem_trace = get_mean_trace(sd,nrem_sel);
    [act_trace,t] = get_mean_trace(sd,act_sel);
end

function [ym,t] = get_mean_trace(sd,sel_ind)

y = {sd.y};
tt = sd(1).t;

mi = min(cellfun(@length,y));
yy = cellfun(@(x) x(1:mi),y,'uni',false);

yy = yy(sel_ind);
yy = [yy{:}]*1000;
t = tt(1:mi)/1000;
ym = median(yy,2);
% Align to trace at 2-3 ms
v = mean(ym(t>=2 & t <=3));
ym = ym-v;