function sdata = get_pooled_ltp_data(keys_pre,keys_post,varargin)
args.suptitle = '';
args.bin_width = 5;
args.rec_gap = 1; % gap between pre and post session
args.slopes_min = 75;
args.slopes_max = 135; % color max
args.tdr_min = 0;
args.tdr_max = 3.5;
% args.exclude_baseline_for_corr = 1;
args.mks = 3; % marker size
args.ch_sel_method = 'ltp_magnitude';
args.plot_corr = 1;
args.corr_cond = {'post'};
args.reverse_col = 0; % reverse it for td ratio and motion, to highlight sleep
args = parseVarArgs(args,varargin{:});
sdata = struct;
% Get the best channel for each pair of pre and post
nExp = length(keys_pre);
kd = struct;
for iExp = 1:nExp
    key_pre = keys_pre(iExp);
    key_post = keys_post(iExp);
    [~,kd(iExp).pre, kd(iExp).post] = get_best_channel(key_pre, key_post,args.ch_sel_method);
end

skeys_pre = [kd.pre];
skeys_post = [kd.post];
sdata.pre = get_pooled_binned_ltpdata(skeys_pre,'slopes',args.bin_width,'invert_time',1);
sdata.post = get_pooled_binned_ltpdata(skeys_post,'slopes',args.bin_width);
% Normalize slopes
sdata = normalize_data_by_pre_mean(sdata);

function dd = normalize_data_by_pre_mean(dd)

n = size(dd.pre.val,2);
mn = nanmean(dd.pre.val,1);
sf = 100./mn;
for i = 1:n
    dd.pre.val(:,i) = dd.pre.val(:,i)*sf(i);
    dd.post.val(:,i) = dd.post.val(:,i)*sf(i);
end
