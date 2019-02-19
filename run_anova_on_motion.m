function run_anova_on_motion(varargin)

args.suptitle = '';
args.bin_width = 5;
args.rec_gap = 1; % gap between pre and post session
args.slopes_min = 75;
args.slopes_max = 135; % color max
args.motion_min = 0;
args.motion_max = 10000;
args.normalize_motion = 1;
args.normalize_avg_motion = 1;
args.tdr_min = 0;
args.tdr_max = 3.5;
% args.exclude_baseline_for_corr = 1;
args.mks = 3; % marker size
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
%% subjects on the rows and time on columns

d = [mdata04.post.val'; mdata06.post.val'];
Shock = [repmat({'04mA'},1,length(keys_post_04)) repmat({'06mA'},1,length(keys_post_06))]';

tab = table(Shock,d(:,1),d(:,2),d(:,3),d(:,4),d(:,5),d(:,6),d(:,7),d(:,8),d(:,9),d(:,10),d(:,11),d(:,12),d(:,13),...
    d(:,14),d(:,15),d(:,16),d(:,17),d(:,18),d(:,19),d(:,20),d(:,21),d(:,22),d(:,23),d(:,24),d(:,25),d(:,26),d(:,27),...
    d(:,28),d(:,29),d(:,30),d(:,31),d(:,32),d(:,33),d(:,34),d(:,35),d(:,36),'VariableNames',{'Shock',...
    't1','t2','t3','t4','t5','t6','t7','t8','t9','t10','t11','t12','t13','t14','t15','t16',...
    't17','t18','t19','t20','t21','t22','t23','t24','t25','t26','t27','t28','t29','t30',...
    't31','t32','t33','t34','t35','t36'});

Time = mdata04.post.t;
rm = fitrm(tab,'t1-t36 ~ Shock','WithinDesign',Time);
ranovatbl = ranova(rm)

%% Pre vs Pre for 0.4 and 0.6mA

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

function dd = normalize_data_by_pre_mean(dd)

n = size(dd.pre.val,2);
mn = nanmean(dd.pre.val,1);
sf = 100./mn;
for i = 1:n
    dd.pre.val(:,i) = dd.pre.val(:,i)*sf(i);
    dd.post.val(:,i) = dd.post.val(:,i)*sf(i);
end