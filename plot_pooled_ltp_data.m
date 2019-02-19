function [h2, sdata, mdata] = plot_pooled_ltp_data(keys_pre,keys_post,keys_day2,varargin)
% function plot_pooled_ltp_data(keys_pre,keys_post,keys_day2,varargin)
% MS 2017-08-12
args.half_errorbars = false;
args.just_get_data_no_plotting = false;
args.suptitle = '';
args.plot_pre = true;
args.plot_tdr = false;
args.remove_nrem = false;
args.remove_sleep = false;
args.keep_just_nrem = false;
args.keep_just_sleep = false;
args.keep_just_motion = false;
args.rec_gap = 0; % gap between pre and post session
args.slopes_min = 75;
args.slopes_max = 135; % color max
args.motion_min = 0;
args.motion_max = 10000;
args.normalize_motion = 1;
args.normalize_avg_motion = 1;
args.tdr_min = 0;
args.tdr_max = 3.5;
args.post_col = [0 0 0];
args.figure = [];
args.slope_axes = [];
% args.exclude_baseline_for_corr = 1;
args.mks = 3; % marker size
args.ch_sel_method = 'ltp_magnitude';
args.plot_corr = 0;
args.plot_just_slope = 0;
args.corr_cond = {'post'};
args.reverse_col = 0; % reverse it for td ratio and motion, to highlight sleep
args.motion_th = 0;
args.motion_th_quantile = 0.5;
args.motion_idx_quantiles = [0 0.95];
args.use_all_data = false;
args.use_pre_event_motion = false;
args.use_pre_event_tdr = false;
args.slope_bw = 5;
args.pre_event_win = 5;
args = parseVarArgs(args,varargin{:});
sdata = struct;
% Get the best channel for each pair of pre and post
nExp = length(keys_pre);
kd = struct;
for iExp = 1:nExp
    key_pre = keys_pre(iExp);
    key_post = keys_post(iExp);
    if ~isempty(keys_day2)
    key_day2 = keys_day2(iExp);
    else
        key_day2 = [];
    end
    [~,kd(iExp).pre, kd(iExp).post, kd(iExp).day2] = get_best_channel(key_pre, key_post,key_day2,args.ch_sel_method,args);
end
fontsz = 12;
skeys_pre = [kd.pre];
skeys_post = [kd.post];
skeys_day2 = [kd.day2];
args.all_motion_th = zeros(1,length(skeys_pre));
nKeys = length(skeys_pre);

for iKey = 1:nKeys
    args.all_motion_idx_bound(iKey,:) = get_common_motion_idx_bound(skeys_pre(iKey),skeys_post(iKey),args);
end

sdata.pre = get_pooled_binned_ltpdata(skeys_pre,'slopes',args,'invert_time',1);
sdata.post = get_pooled_binned_ltpdata(skeys_post,'slopes',args,'invert_time',0);
sdata.day2 = get_pooled_binned_ltpdata(skeys_day2,'slopes',args,'invert_time',0);

% Normalize slopes
sdata = normalize_data_by_pre_mean(sdata);



% Motion
mdata.pre = get_pooled_binned_ltpdata(skeys_pre,'motion',args,'invert_time',1);
mdata.post = get_pooled_binned_ltpdata(skeys_post,'motion',args,'invert_time',0);
mdata.day2 = get_pooled_binned_ltpdata(skeys_day2,'motion',args,'invert_time',0);
% Normalize if asked for
if args.normalize_motion
    mdata = normalize_data_by_pre_mean(mdata);
end

if args.just_get_data_no_plotting
    h2 = [];
    return
end
% Theta delta ratio
tddata.pre = get_pooled_binned_ltpdata(skeys_pre,'tdratio',args,'invert_time',1);
tddata.post = get_pooled_binned_ltpdata(skeys_post,'tdratio',args,'invert_time',0);
tddata.day2 = get_pooled_binned_ltpdata(skeys_day2,'tdratio',args,'invert_time',0);


% Finally, time to plot things
args.main_data_p_col = 5;

if args.plot_just_slope
    gs = [2,args.main_data_p_col];
else
    gs = [6,args.main_data_p_col];
    if ~args.plot_tdr
        gs = [4,args.main_data_p_col];
    end
end
if ~isempty(keys_day2)
    gs(2) = gs(2)+2;
end
% if isempty(args.figure)
%     fig1 = 100;
% else
%     fig1 = args.figure;
% end
% figure(fig1)
figure
if args.plot_just_slope
    set(gcf,'Position',[344,335,1023,633],'color','w')
else
   set(gcf,'Position',[2104,336,961,1021],'color','w')
end
dcol = 1:args.main_data_p_col;

%% Slopes
% Raw normalized slopes of all mice

h1 = msubplot(1,dcol,gs);
rev_col = false;
if ~isempty(keys_day2)
    colbar = 0;
else
    colbar = 1;
end
plot_raw_data(sdata,h1,args.slopes_min,args.slopes_max,rev_col,colbar)
title('Norm epsp slope')
ylabel('Mouse #')

if ~isempty(keys_day2)
    hr2 = msubplot(1,args.main_data_p_col+(1:2),gs);

plot_raw_data_day2(sdata,hr2,args.slopes_min,args.slopes_max,rev_col)
% title('Day2')
end
% Day 2 data


if args.plot_just_slope
set(gca,'Fontsize',fontsz)
end
if isempty(args.slope_axes)
h2 = msubplot(2,dcol,gs);
else
    h2 = args.slope_axes;
end
% Averaged across mice
plot_avg_data(sdata,h2,args.slopes_min,args.slopes_max,args)
ylabel('Avg norm slope')
xlabel('Time (min)')
if args.plot_just_slope
set(gca,'Fontsize',fontsz)
end
if ~isempty(keys_day2)
    hd2 = msubplot(2,args.main_data_p_col+(1:2),gs);

plot_avg_data_day2(sdata,hd2,args.slopes_min,args.slopes_max,args)
end
% Day 2 data
if ~args.plot_just_slope
    %% Motion
    h3 = msubplot(3,dcol,gs);
    rev_col = false;
    mmax = quantile([mdata.pre.val(:); mdata.post.val(:)],0.975);
    plot_raw_data(mdata,h3,args.motion_min,mmax,rev_col,1)
    ylabel('Mouse #')
    title('Motion index')
    
    h4 = msubplot(4,dcol,gs);
    % Averaged across mice
    if args.normalize_motion || args.normalize_avg_motion
        ylab = 'Norm Avg Motion Index';
        vmin = 0;
        vmax = 210;
    else
        ylab = 'Avg motion index';
        vmin = args.motion_min;
        vmax = args.motion_max;
    end
    
    if args.normalize_avg_motion
        plot_norm_avg_data(mdata,h4,vmin,vmax,args)
    else
        plot_avg_data(mdata,h4,vmin,vmax,args)
    end
    ylabel(ylab)
    
    
    
    %% td-ratio
    if args.plot_tdr
    h5 = msubplot(5,dcol,gs);
    rev_col = false;
    tdmax = quantile([tddata.pre.val(:); tddata.post.val(:)],0.975);
    plot_raw_data(tddata,h5,args.tdr_min,tdmax,rev_col,1)
    ylabel('Mouse #')
    title('\theta/\delta ratio')
    
    h6 = msubplot(6,dcol,gs);
    % Averaged across mice
    plot_avg_data(tddata,h6,args.tdr_min,args.tdr_max,args)
    ylabel('Avg \theta/\delta ratio')
    xlabel('Time (min)')
    
    % Plot all correlation data now
%     plot_corr_data(sdata,mdata,tddata,args.corr_cond,args)
%     figure(fig1)
    end
end
ms_suptitle(args.suptitle,'yPosition',0.975)


function plot_norm_avg_data(dd,h,vmin,vmax,args)
col = 'k';
ad1 = dd.pre.val;
ad2 = dd.post.val;
ad3 = dd.day2.val;
t1 = dd.pre.t;
t2 = dd.post.t;
t3 = dd.day2.t;

if args.rec_gap > 0
    t2 = t2 + args.rec_gap;
end
t3 = t3+t2(end);

nMice = size(ad1,2);

% average now
mad1 = nanmean(ad1,2);
se1 = nanstd(ad1,[],2)/sqrt(nMice);
mad2 = nanmean(ad2,2);
se2 = nanstd(ad2,[],2)/sqrt(nMice);
mad3 = nanmean(ad3,2);
se3 = nanstd(ad3,[],2)/sqrt(nMice);

% Normalize the average now
amad1 = mean(mad1);
sf = 100/amad1;
mad1 = mad1 * sf;
se1 = se1 * sf;
mad2 = mad2 * sf;
se2 = se2 * sf;
mad3 = mad3 * sf;
se3 = se3 * sf;

axes(h)
ebneg1 = se1;
ebneg2 = se2;
ebneg3 = se3;
if args.half_errorbars
    ebneg1 = zeros(size(se1));
    ebneg2 = zeros(size(se1));
    ebneg3 = zeros(size(se3));
end
plot(t1,mad1,'O','color',col,'markersize',args.mks,'markerfacecolor',col)
hold on
errorbar(t1,mad1,ebneg1,se1,'color',col,'linestyle','none')
plot(t2,mad2,'O','color',col,'markersize',args.mks,'markerfacecolor',col)
errorbar(t2,mad2,ebneg2,se2,'color',col,'linestyle','none')
plot(t3,mad3,'O','color',col,'markersize',args.mks,'markerfacecolor',col)
errorbar(t3,mad3,ebneg3,se3,'color',col,'linestyle','none')


baseline = mean(mad1);
plot(xlim,[baseline baseline],'r','linewidth',1)
plot([0 0],ylim,'color',[0.8 0.8 0.8])
axis tight
xlim([t1(1) t2(end)])
ylim([vmin vmax])
box off
grid on
grid minor

function plot_corr_data(sdata,mdata,tddata,corr_cond,args)

% Plot pre and post separately
% Motion versus LTP - pool pre and post
stypes = {'pre','post'};
ttypes = {'Before','After'};
nCond = length(corr_cond);
% off = args.main_data_p_col;
figure
set(gcf,'Position',[2818,254,492,1098],'color','w')
gs = [3 1];
for iCond = 1:nCond
    %     pcol = off + 1 +((iCond-1):(iCond));
    %     pcol = off+iCond;
    pcol = iCond;
    % Motion versus LTP
    h = msubplot(1,pcol,gs);
    ccond = corr_cond{iCond};
    X1 = mdata.(ccond).val(:);
    Y = sdata.(ccond).val(:);
    if args.normalize_motion
        mlab = 'Norm Motion Index';
    else
        mlab = 'Motion Index';
    end
    ylab = 'Norm epsp slope';
    titstr = ttypes{strcmp(ccond,stypes)};
    plot_corr(X1,Y,h,mlab,ylab,titstr);
    
    
    %% TDratio versus LTP
    h = msubplot(2,pcol,gs);
    X2 = tddata.(ccond).val(:);
    % TD ratio versus slopes
    xlab = '\theta/\delta ratio';
    ylab = 'Norm epsp slope';
    plot_corr(X2,Y,h,xlab,ylab);
    
    %% Plot motion versus th/del ratio
    h = msubplot(3,pcol,gs);
    ylab = '\theta/\delta ratio';
    plot_corr(X1,X2,h,mlab,ylab)
end
ms_suptitle(args.suptitle,'yPosition',0.975)
%%
function plot_corr(X,Y,h,xlab,ylab,titstr)
if nargin < 6
    titstr = '';
end
axes(h)
lm = fitlm(X,Y,'linear','RobustOpts','on');
y = predict(lm,X);
p = coefTest(lm);
rs = roundn(lm.Rsquared.Ordinary,-2);
plot(X,Y,'k.','markerfacecolor','none')
hold on
plot(X,y,'r','linewidth',2)
axis tight square
grid on; box off
title(sprintf('%s\n%s %s',titstr,get_plessthan_str(p),['R^{2} = ' num2str(rs)]),'fontweight','normal','fontsize',10)
xlabel(xlab)
ylabel(ylab)
xlim(quantile(X(:),[0 0.975]))
ylim(quantile(Y(:),[0 0.975]))
%%
function plot_avg_data_day2(dd,h,vmin,vmax,args)

precol = 'k';

if ~args.plot_pre
    dd.pre.val = nan(size(dd.pre.val));
end
ad1 = dd.pre.val;
ad3 = dd.day2.val;
t3 = dd.day2.t;

nMice = size(ad3,2);

% average now
% average now
mad1 = nanmean(ad1,2);
mad3 = nanmean(ad3,2);
se3 = nanstd(ad3,[],2)/sqrt(nMice);


axes(h)

plot(t3,mad3,'O','color',precol,'markersize',args.mks,'markerfacecolor',precol)
hold on
errorbar(t3,mad3,se3,'color',precol,'linestyle','none')


baseline = mean(mad1);
plot(xlim,[baseline baseline],'r','linewidth',1)
plot([0 0],ylim,'color',[0.8 0.8 0.8])
xtstr = [0 50 100];
xtloc = xtstr;
set(gca,'XTick',xtloc,'XTickLabel',xtstr,'FontSize',10,'yticklabel',[])
axis tight
xlim([0 t3(end)])
ylim([vmin vmax])
box off
grid on
grid minor
%%
function plot_avg_data(dd,h,vmin,vmax,args)
precol = 'k';
postcol = args.post_col;
if ~args.plot_pre
    dd.pre.val = nan(size(dd.pre.val));
end
ad1 = dd.pre.val;
ad2 = dd.post.val;
% ad3 = dd.day2.val;
t1 = dd.pre.t;
t2 = dd.post.t;
% t3 = dd.day2.t;
if args.rec_gap > 0
    t2 = t2 + args.rec_gap;
end
% t3 = t3 + t2(end);
nMice = size(ad1,2);

% average now
mad1 = nanmean(ad1,2);
se1 = nanstd(ad1,[],2)/sqrt(nMice);
mad2 = nanmean(ad2,2);
se2 = nanstd(ad2,[],2)/sqrt(nMice);
% mad3 = nanmean(ad3,2);
% se3 = nanstd(ad3,[],2)/sqrt(nMice);


axes(h)
plot(t1,mad1,'O','color',precol,'markersize',args.mks,'markerfacecolor',precol)
hold on
ebneg1 = se1;
ebneg2 = se2;
if args.half_errorbars
    ebneg1 = zeros(size(se1));
    ebneg2 = zeros(size(se2));
end
errorbar(t1,mad1,ebneg1,se1,'color',precol,'linestyle','none')
tadj = t2+args.rec_gap;
plot(tadj,mad2,'O','color',postcol,'markersize',args.mks,'markerfacecolor',postcol)
errorbar(tadj,mad2,ebneg2,se2,'color',postcol,'linestyle','none')

% plot(t3,mad3,'O','color',postcol,'markersize',args.mks,'markerfacecolor',postcol)
% errorbar(t3,mad3,se3,'color',postcol,'linestyle','none')


baseline = mean(mad1);
plot(xlim,[baseline baseline],'r','linewidth',1)
plot([0 0],ylim,'color',[0.8 0.8 0.8])
xtstr = [-100 -50 0 50 100 150];
xtloc = [-100 -50 0 [50 100 150]+args.rec_gap];
set(gca,'XTick',xtloc,'XTickLabel',xtstr,'FontSize',10)
axis tight
xlim([t1(1) t2(end)])
ylim([vmin vmax])
box off
grid on
grid minor
%%
function plot_raw_data(dd,h,cmin,cmax,rev_col,colbar)
axes(h)
% Stitch pre and post
C = [dd.pre.val; dd.post.val]';
x = 1:size(C,1);
t = [dd.pre.t; dd.post.t]';
colormap(h,'default')
if rev_col
    colormap(h,flipud(colormap))
end
% pcolor(x,y,C)
h = imagesc(t,x,C);
set(h,'alphadata',~isnan(C))
set(gca, 'CLim', [cmin, cmax],'FontSize',10,'ytick',x)
if colbar
colorbar
end


function plot_raw_data_day2(dd,h,cmin,cmax,rev_col)
axes(h)
% Stitch pre and post
C = [dd.day2.val]';
x = 1:size(C,1);
t = [dd.day2.t]';
colormap(h,'default')
if rev_col
    colormap(h,flipud(colormap))
end
% pcolor(x,y,C)
h = imagesc(t,x,C);
set(h,'alphadata',~isnan(C))
set(gca, 'CLim', [cmin, cmax]);
colorbar
set(gca,'yticklabel',[],'FontSize',10)


function dd = normalize_data_by_pre_mean(dd)

n = size(dd.pre.val,2);
mn = nanmean(dd.pre.val,1);
sf = 100./mn;
for i = 1:n
    dd.pre.val(:,i) = dd.pre.val(:,i)*sf(i);
    dd.post.val(:,i) = dd.post.val(:,i)*sf(i);
    if ~isempty(dd.day2.val)
    dd.day2.val(:,i) = dd.day2.val(:,i)*sf(i);
    else
        dd.day2.val = [];
    end
end



