function [fkeys_pre, fkeys_post, p] = filter_keys_by_stability(keys_pre, keys_post)
% filter_keys_by_stability(keys_pre, keys_post)
% Use bootstrapping and remove sessions that didn't have a 'stable'
% baseline during the 'pre' condition.
% Mani Subramaniyan 2019-02-19
%
alpha = 0.05;
tmp = struct;
[~,tmp.sdata] = plot_pooled_ltp_data(keys_pre,keys_post,[],'slope_bw',1,'use_pre_event_motion',2,'motion_th_quantile',0.5,...
    'plot_just_slope',1,'motion_idx_quantiles',[0 0.99],'use_pre_event_tdr',0,'pre_event_win',2,'use_all_data',0,'just_get_data_no_plotting',1);
%% get p value
for i = 1:size(tmp.sdata.pre.val,2)
    tmp.p(i) = test_baseline_stability_kolmogorov_dimitri(tmp.sdata.pre.val(:,i),'debug',0);
end
p = tmp.p;
%%
fkeys_pre = keys_pre(tmp.p > alpha);
fkeys_post = keys_post(tmp.p > alpha);
