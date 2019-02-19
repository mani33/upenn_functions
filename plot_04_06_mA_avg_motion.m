function plot_04_06_mA_avg_motion(varargin)

args.suptitle = '';
args.bin_width = 5;
args.rec_gap = 5; % gap between pre and post session
args.slopes_min = 75;
args.slopes_max = 135; % color max
args.motion_min = 0;
args.motion_max = 10000;
args.normalize_motion = 1;
args.normalize_avg_motion = 1;
args.tdr_min = 0;
args.tdr_max = 3.5;
% args.exclude_baseline_for_corr = 1;
args.mks = 5; % marker size
args.ch_sel_method = 'ltp_magnitude';
args.plot_corr = 1;
args.corr_cond = {'post'};
args.reverse_col = 0; % reverse it for td ratio and motion, to highlight sleep
args = parseVarArgs(args,varargin{:});

datatype = 'motion';
fn = 'z:\users\mani\matlab\analysis\FearConditioning\Sessions for Analysis_2017_07_22.xlsx';
% 0.4 mA shock
keys_pre_04 = get_keys_from_xlsrange(fn,1,'G2:G10'); 
keys_post_04 = get_keys_from_xlsrange(fn,1,'H2:H10');

% 0.6 mA shock
keys_pre_06 = get_keys_from_xlsrange(fn,1,'G11:G19'); 
keys_post_06 = get_keys_from_xlsrange(fn,1,'H11:H19');

[keys_pre_04, keys_post_04] = filter_best_keys(keys_pre_04,keys_post_04,args);
[keys_pre_06, keys_post_06] = filter_best_keys(keys_pre_06,keys_post_06,args);

% 0.4 versus 0.6 mA during post session
mdata04.pre = get_pooled_binned_ltpdata(keys_pre_04,datatype,args.bin_width,'invert_time',1);
mdata04.post = get_pooled_binned_ltpdata(keys_post_04,datatype,args.bin_width);

mdata06.pre = get_pooled_binned_ltpdata(keys_pre_06,datatype,args.bin_width,'invert_time',1);
mdata06.post = get_pooled_binned_ltpdata(keys_post_06,datatype,args.bin_width);

if strcmp(datatype,'slopes')
    mdata04 = normalize_data_by_pre_mean(mdata04);
    mdata06 = normalize_data_by_pre_mean(mdata06);
end




figure
h = subplot(1,1,1);
ph = nan(1,2);
ph(1) = plot_norm_avg_data(mdata04,h,0, 200,args,'k');
ph(2) = plot_norm_avg_data(mdata06,h,0, 200,args,'r');
leg = legend(ph,{sprintf('0.4 mA (n = %u)',length(keys_pre_04)),sprintf('0.6 mA (n = %u)',length(keys_pre_06))});
set(leg,'box','off')
xlabel('Time (min)')
ylabel('Avg norm motion index')

function plot_avg_data(dd,h,vmin,vmax,args,col)

ad1 = dd.pre.val;
ad2 = dd.post.val;
t1 = dd.pre.t;
t2 = dd.post.t;

if args.rec_gap > 0
    t2 = t2 + args.rec_gap;
end

nMice = size(ad1,2);

% average now
mad1 = nanmean(ad1,2);
se1 = nanstd(ad1,[],2)/sqrt(nMice);
mad2 = nanmean(ad2,2);
se2 = nanstd(ad2,[],2)/sqrt(nMice);

axes(h)
plot(t1,mad1,'O','color',col,'markersize',args.mks,'markerfacecolor',col)
hold on
errorbar(t1,mad1,se1,'color',col,'linestyle','none')
plot(t2,mad2,'O','color',col,'markersize',args.mks,'markerfacecolor',col)
errorbar(t2,mad2,se2,'color',col,'linestyle','none')
baseline = mean(mad1);
plot(xlim,[baseline baseline],'r','linewidth',1)
plot([0 0],ylim,'color',[0.8 0.8 0.8])
axis tight
xlim([t1(1) t2(end)])
ylim([vmin vmax])
box off
grid on
grid minor


function [skeys_pre, skeys_post] = filter_best_keys(keys_pre,keys_post,args)
nExp = length(keys_pre);
kd = struct;
for iExp = 1:nExp
    key_pre = keys_pre(iExp);
    key_post = keys_post(iExp);
    [~,kd(iExp).pre, kd(iExp).post] = get_best_channel(key_pre, key_post,args.ch_sel_method);
end

skeys_pre = [kd.pre];
skeys_post = [kd.post];


function ph = plot_norm_avg_data(dd,h,vmin,vmax,args,col)

ad1 = dd.pre.val;
ad2 = dd.post.val;
t1 = dd.pre.t;
t2 = dd.post.t;

if args.rec_gap > 0
    t2 = t2 + args.rec_gap;
end

nMice = size(ad1,2);

% average now
mad1 = nanmean(ad1,2);
se1 = nanstd(ad1,[],2)/sqrt(nMice);
mad2 = nanmean(ad2,2);
se2 = nanstd(ad2,[],2)/sqrt(nMice);
% Normalize the average now
amad1 = mean(mad1);
sf = 100/amad1;
mad1 = mad1 * sf;
se1 = se1 * sf;
mad2 = mad2 * sf;
se2 = se2 * sf;

axes(h)
ph = plot(t1,mad1,'O-','color',col,'markersize',args.mks,'markerfacecolor',col);
hold on
errorbar(t1,mad1,se1,'color',col,'linestyle','none')
plot(t2,mad2,'O-','color',col,'markersize',args.mks,'markerfacecolor',col)
errorbar(t2,mad2,se2,'color',col,'linestyle','none')
baseline = mean(mad1);
plot([-120 180],[baseline baseline],'b--','linewidth',1)
plot([0 0],ylim,'color',[0.8 0.8 0.8])
axis tight
xlim([t1(1) t2(end)])
ylim([vmin vmax])
box off
grid on
grid minor
