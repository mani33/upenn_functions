function plot_avg_ltp(keys_pre,keys_post,h,varargin)
% plot ltp data - pre and post averaged over all keys (keys must be for
% individual selected channels, one for each mouse.
args.mouse141key2 = [];
args.plot_motion = true;
args = parseVarArgs(args,varargin{:});

if isempty(h)
    h = gca;
end
yl = [75 150];
bin_width = 5;
npre = length(keys_pre);
npost = length(keys_post);
assert(npre==npost,'supply equal number of pre and post keys')
n = npre;
tmp = struct;
pc = 'k';
ms = 6;
% Compute average pre
for i = 1:n
    key_pre = keys_pre(i);
    [t_minutes,y] = get_slope_data(key_pre);
    [resp,~,tb] = bin_slopes(y,t_minutes,bin_width);
    tmp.sf(i) = 100/mean(resp);
    nresp = resp*tmp.sf(i);
    tmp.ori_pre{i} = nresp(:);
end
% Pad nan's for sessions that had 45 min baseline instead of 2h
% mpulses = (120/bin_width);
mpulses = median(cellfun(@length,tmp.ori_pre));
for i = 1:n
    slopes = tmp.ori_pre{i};
    np = length(slopes);
    if np < mpulses
        tmp.slopes_pre(:,i) = [nan(mpulses-np,1); slopes];
    else
        tmp.slopes_pre(:,i) = slopes;
    end
end
mean_pre = nanmean(tmp.slopes_pre,2);
% Compute std error
s = nanstd(tmp.slopes_pre,[],2);
se = s/sqrt(n);
% axes(h)
subplot(3,1,1)
plot([-122 182],[100 100],'r-','linewidth',1)
hold on
ttpre = tb-tb(end);
plot(ttpre,mean_pre,'O','color',pc,'markersize',ms,'markerfacecolor',pc)
hold on
errorbar(ttpre,mean_pre,se,'color',pc,'linestyle','none')


% post
for i = 1:n
    key_post = keys_post(i);
    if key_post.animal_id == 141 && ~isempty(args.mouse141key2)
        [t_minutes,y] = get_slope_data_141(key_post,args.mouse141key2);
    else
        [t_minutes,y] = get_slope_data(key_post);
    end
    [resp,~,tt] = bin_slopes(y,t_minutes,bin_width);
    tmp.slopes_post(:,i) = resp*tmp.sf(i);
end
mean_post = nanmean(tmp.slopes_post,2);

% Compute std error
s = nanstd(tmp.slopes_post,[],2);
se = s/sqrt(n);
hold on
plot(tt,mean_post,'O','color',pc,'markersize',ms,'markerfacecolor',pc)
hold on
errorbar(tt,mean_post,se,'color',pc,'linestyle','none')
set(gca,'FontSize',16)
plot([0 0],ylim,'k--')
xlim([ttpre(1)-2 inf])
ylim(yl)
plot([0 0],ylim,'k--')

% PLot motion
[m_data,tm,td_data,t_td] = get_mouse_motion_and_tdratio_for_ltp(keys_pre,keys_post,'avg',true);
if args.plot_motion
    subplot(3,1,2)
    plot(tm,m_data,'k.')
    axis tight
    subplot(3,1,3)
    plot(t_td,td_data,'k.')
    hold on
    axis tight
    plot(xlim,[1 1],'r','linewidth',2)
    
end
% Make a color plot
figure
prepost = [tmp.slopes_pre(end-8:end,:); tmp.slopes_post]';
tprepost = [ttpre(end-8:end) tt];
imagesc(tprepost,1:n,prepost)
h = colorbar;
set(gca,'FontSize',16)
hold on
plot([0 0],ylim,'w--','linewidth',2)
xlabel('Time')
ylabel('Mouse #')


% correlation between ltp and freezing
post_ltp = nanmean(tmp.slopes_post,1);
figure; 
lowcol = [0.6 0.6 0.6];
for i = 1:n
    try
        tmp.freezing(i) = fetch1(beh.FearCond(keys_pre(i)),'percent_freezing');
        tmp.shock(i) = fetch1(beh.FearCond(keys_pre(i)),'shock_intensity');
    catch
        tmp.freezing(i)   = nan;
        tmp.shock(i) = nan;
    end
if tmp.shock(i)==0.4
    col = lowcol;
else
    col = 'k';
end
plot(tmp.freezing(i),post_ltp(i),'O','color',col,'markerfacecolor',col,'markersize',12)
hold on
end
xlim([0 100])
ylim([90 145])
grid on
set(gca,'FontSize',16)
hold on
plot([0 15],[100 100],'k--','linewidth',2)
xlabel('Freezing (%)')
ylabel('Mean LTP')
text(10,140,'0.4 mA','color',lowcol,'fontsize',16,'fontweight','bold')
text(10,135,'0.6 mA','color','k','fontsize',16,'fontweight','bold')
box off




function [t_minutes,y] = get_slope_data_141(key_post1,key_post2)

% Post-shock
[t,slopes1] = fetchn(cstim.FepspSlope(key_post1),'event_ts','fepsp_slope');

sel = 1:100;
t_minutes = format_time(t);
t_1 = t_minutes(sel);
slopes1 = slopes1(sel);

% Get second session for stitching
[t,slopes2] = fetchn(cstim.FepspSlope(key_post2),'event_ts','fepsp_slope');
t_2 = format_time(t)+120;

t_minutes = [t_1(:); (100:119)';t_2(:)];
y = [slopes1(:); nan(20,1); slopes2(:)];

function [t_minutes,slopes] = get_slope_data(key)
[t,slopes] = fetchn(cstim.FepspSlope(key),'event_ts','fepsp_slope');

sel = 1:min([180 length(slopes)]);
if length(sel)~=length(slopes)
    warning('Number of pulses were %u',length(slopes))
end
t_minutes = format_time(t);
t_minutes = t_minutes(sel);
slopes = slopes(sel);