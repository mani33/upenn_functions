function motion_idx_bound = get_common_motion_idx_bound(ckey_pre,ckey_post,args)
% motion_idx_bound = get_common_motion_idx_bound(ckey_pre,ckey_post,args)
% Mani 2018-05-04
% Get common motion index boundary for before and after shock conditions.
if args.use_all_data
    motion_idx_bound = [-inf inf];
    return
end
q = args.motion_th_quantile;
winstr = sprintf('win = %u',args.pre_event_win);
mdx_pre = fetchn(cont.PreEventMotion(ckey_pre,winstr),'dist_var');
if isempty(mdx_pre)
    mdx_pre = fetchn(cont.PreEventMotionInVid(ckey_pre,winstr),'dist_var');
end
if ~isempty(ckey_post)
mdx_post = fetchn(cont.PreEventMotion(ckey_post,winstr),'dist_var');
if isempty(mdx_post)
    mdx_post = fetchn(cont.PreEventMotionInVid(ckey_post,winstr),'dist_var');
end
else
    mdx_post = NaN;
end
qpre = quantile(mdx_pre,q);
qpost = quantile(mdx_post,q);

qthresh = min([qpre qpost]);


% Pick only those motion indices that are above the bigger of the pre vs
% post threshold
mpre = mdx_pre(mdx_pre >= qthresh);
mpost = mdx_post(mdx_post >= qthresh);

% From these selected population, find a common range that accommodates
% most of the points for both pre and post.
qa = args.motion_idx_quantiles;
qqpre = quantile(mpre,qa);
qqpost = quantile(mpost,qa);

motion_idx_bound = [max([qqpre(1) qqpost(1)]), min([qqpre(2) qqpost(2)])];
fprintf('Selecting common motion index bound for pre and post session for mouse: %u\n',ckey_pre.animal_id)
