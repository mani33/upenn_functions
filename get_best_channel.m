function [ch_ind, best_pre_key, best_post_key, best_day2_key]  = get_best_channel(key_pre,key_post,key_nextday_baseline,judge,args)
% [ch_ind, best_pre_key, best_post_key]  = get_best_channel(key_pre,key_post,key_nextday_baseline,judge,varargin)
% MS 2017-08-12
% MS 2018-01-18

winstr = sprintf('win = %u',args.pre_event_win);
pckeys_pre = fetch(cstim.SlopeEg & cstim.FepspSlope(key_pre,'smooth_method_num = 0'));
pckeys_post = fetch(cstim.SlopeEg & cstim.FepspSlope(key_post,'smooth_method_num = 0'));
if ~isempty(key_nextday_baseline)
    pckeys_nextday = fetch(cstim.SlopeEg & cstim.FepspSlope(key_nextday_baseline,'smooth_method_num = 0'));
else
    pckeys_nextday = [];
end
nKeys = length(pckeys_pre);

% Sort keys based on channel number so they match for pre and post
% condition
chnpre = [pckeys_pre.chan_num];
chnpost = [pckeys_post.chan_num];
if ~isempty(pckeys_nextday)
    chnday2 = [pckeys_nextday.chan_num];
else
    chnday2 = [];
end
common = intersect(chnpre,chnpost);
for iC = 1:length(common)
    ckeys_pre(iC) = pckeys_pre(chnpre==common(iC));
    ckeys_post(iC) = pckeys_post(chnpost==common(iC));
    if ~isempty(chnday2)
        ckeys_nextday(iC) = pckeys_nextday(chnday2==common(iC));
    else
        ckeys_nextday = [];
    end
end


if nKeys == 1
    ch_ind = 1;
    best_pre_key = ckeys_pre;
    best_post_key = ckeys_post;
    best_day2_key = ckeys_nextday;
    return
end
ltp = nan(1,nKeys);

for i = 1:nKeys
    kpre = ckeys_pre(i);
    kpost = ckeys_post(i);
    %
    fprintf('Selecting common motion index bound for pre and post session for key: %u\n',i)
    mib = get_common_motion_idx_bound(kpre,kpost,args);
    
    mdx_pre = fetchn(cont.PreEventMotion(kpre,winstr),'dist_var');
    
    if isempty(mdx_pre)
        mdx_pre = fetchn(cont.PreEventMotionInVid(kpre,winstr),'dist_var');
    end
    mdx_post = fetchn(cont.PreEventMotion(kpost,winstr),'dist_var');
    if isempty(mdx_post)
        mdx_post = fetchn(cont.PreEventMotionInVid(kpost,winstr),'dist_var');
    end
    [slopes_pre,tpre] = fetchn(cstim.FepspSlope(kpre,'smooth_method_num = 0'),'fepsp_slope','event_ts');
    sel = mdx_pre > mib(1) & mdx_pre < mib(2);
    slopes_pre(~sel) = nan;
    
    [slopes_post] = fetchn(cstim.FepspSlope(kpost,'smooth_method_num = 0'),'fepsp_slope');
    sel = mdx_post > mib(1) & mdx_post < mib(2);
    slopes_post(~sel) = nan;
    
%     tpre = double((tpre-tpre(1)))*(1e-6)/60;
    val = slopes_post/nanmean(slopes_pre);
    ltp(i) = nanmean(val);
%     if args.use_baseline_stability
%         p = test_for_baseline_stability(tpre,slopes_pre);
%         % Higher the p value, the more insignificant is the slope value of
%         % the baseline data points. So, we give some weight to that.
%         ltp(i) = ltp(i) + p;
%     end
    
end


switch judge
    case 'ltp_magnitude'
        measure = ltp;
    otherwise
        error('illegal judging method supplied')
end

[~,ch_ind] = max(measure);

best_pre_key = ckeys_pre(ch_ind);
best_post_key = ckeys_post(ch_ind);
if ~isempty(key_nextday_baseline)
    best_day2_key = ckeys_nextday(ch_ind);
else
    best_day2_key = [];
end
if best_pre_key.session_start_time == 5336101341
    best_pre_key.chan_num = 14;
    best_post_key.chan_num = 14;
    % We do this because the last 5 points or so in the baseline ended up
    % unusually flat. this changes the overall nanmean to weird location.
end

fprintf('Best channel for mouse %u is %u\n',best_pre_key.animal_id,best_pre_key.chan_num)

