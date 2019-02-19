function plot_pre_vs_post_motion(keys_pre,keys_post,varargin)
% function plot_pre_vs_post_motion(key_pre,key_post,varargin)
% Mani 2018-05-05

args.motion_th_quantile = 0.5;
args.motion_idx_quantiles = [0 0.99];

args = parseVarArgs(args,varargin{:});
nKeys = length(keys_pre);
d = struct;
for iKey = 1:nKeys
    key_pre = keys_pre(iKey);
    key_post = keys_post(iKey);
    mib = get_common_motion_idx_bound(key_pre,key_post,args);
    
    mdx_pre = fetchn(cont.PreEventMotion(key_pre),'dist_var');
    if isempty(mdx_pre)
        mdx_pre = fetchn(cont.PreEventMotionInVid(key_pre),'dist_var');
    end
    mdx_post = fetchn(cont.PreEventMotion(key_post),'dist_var');
    if isempty(mdx_post)
        mdx_post = fetchn(cont.PreEventMotionInVid(key_post),'dist_var');
    end
   
    mpre = mdx_pre(mdx_pre >= mib(1) & mdx_pre < mib(2));
    mpost = mdx_post(mdx_post >= mib(1) & mdx_post < mib(2));    
    d.mean(iKey,1:2) = [nanmean(mpre)  nanmean(mpost)];
    d.median(iKey,1:2) = [nanmedian(mpre) nanmedian(mpost)];
end
clf
subplot(1,2,1)
title('Mean')
plot_bar_a_vs_b(d.mean(:,1),d.mean(:,2),{'pre','post'},'show_legend',0,'paired',1,'boxplot',0)

subplot(1,2,2)
title('Median')
plot_bar_a_vs_b(d.median(:,1),d.median(:,2),{'pre','post'},'show_legend',0,'paired',1,'boxplot',0)
