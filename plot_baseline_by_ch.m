function plot_baseline_by_ch(key,varargin)
% Plot fepsp response measure along with motion, theta delta ratio


args.plot_his = 1;
args.minfo_exists = 1; % motion info exists or not
args.bin_motion = 1;
args.plot_avg_trace = 1;
args.plot_motion = 1;
args.plot_tdratio = 1;
args.plot_corr = 1;
args.overlay_motion = 1;
args.overlay_tdratio = 1;
args.max_ch = inf;
args.bw = 5; % in minutes for binning slopes and motion and tdratio
args = parseVarArgs(args,varargin{:});

ch_keys = fetch(cstim.SlopeEg & cstim.FepspSlope(key),'smooth_method_num = 0');

nKeys = min(length(ch_keys),args.max_ch);


havg = nan(1,2);
m = 1;
if args.plot_his
    m = 2;
end

if args.plot_avg_trace
    m = 3;
end
% xmax = 180;

pc = 'k';
plot_raw_resp = 1;
plot_binned = 1;
yl = [25 175];
coal = 2;
gs = [nKeys+2 m+coal];
h = nan(1,nKeys);
for i = 1:nKeys
    h(i) = msubplot(i,1:(coal+1),gs);
    % Baseline
    ch_key = ch_keys(i);
    if i == 1
        tit_str = sprintf('Mouse # %u: ch %u', ch_key.animal_id,ch_key.chan_num);
    else
        tit_str = sprintf('ch %u', ch_key.chan_num);
    end
    
    [t,binnedSlopes] = fetchn(cstim.SlopeBinned(ch_key,'smooth_method_num = 0',sprintf('slope_bw = %0.2f',args.bw)),'event_ts','fepsp_slope');
    sfac = 100/mean(binnedSlopes);
    binnedSlopes = binnedSlopes * sfac;
    t_minutes = format_time(t);
    ttt = t_minutes;
    %     ttt = t_minutes-max(t_minutes);
    plot_epsp_resp(ttt,binnedSlopes,plot_raw_resp,plot_binned,bw,h(i),pc,[],[],[])
    set(gca,'FontSize',12)
    [slopes_pre,~,~] = bin_mouse_data(binnedSlopes,ttt,bw);
    
    plot([0 0],ylim,'k--')
    %     xl1 = min([-46 ttt(1)]);
    %     xlim([xl1 xmax])
    grid on
    cylabel('Norm fepsp slope',i==1)
    
    % Plot avg trace
    h2 = msubplot(i,4,gs);
    havg(1) = plot_avg_trace(ch_key,h2,'k');
    xlabel('Time (ms)')
    title(tit_str)
    % Plot histology
    if args.plot_his
        h3 = msubplot(i,5,gs);
        plot_histology(ch_key,1,'r',h3)
    end
end

if args.minfo_exists
    ms = 4;
    [t, m_data] = get_motion_data(ch_key);
    tm = t-t(1);
    [bm_data,~,btm] = bin_mouse_data(m_data,tm,bw);
end
[t, td_data] = get_tdratio_data(ch_key);
t_td = t - t(1);

% [m_data,tm,td_data,t_td] = get_mouse_motion_and_tdratio_for_ltp(key_pre,key_post);
[btd_data,~,bt_td] = bin_mouse_data(td_data,t_td,bw);



% Plot motion data

datarange = 100;
mth = nan;
mstr = '';
if args.minfo_exists
    hm = [];
    if args.plot_motion
        ms = 2;
        i = i+1;
        hm = msubplot(i,1:3,gs);
        plot(tm,m_data,'O','color',pc,'markersize',ms,'markerfacecolor',pc)
        hold on
        plot([0 0],ylim,'k--')
        grid on
        linkaxes([h hm],'x')
        %     xlim([xl1 xmax])
        ylim([0 quantile(m_data,0.99)])
        ylabel('Motion Index')
        box off
        
        
        
        
        % Overlay the motion data on the neural responses
       
        if args.overlay_motion
            rm = bm_data-nanmin(bm_data);
            rm = datarange*rm/nanmax(rm);
            rm = 100+(rm-nanmean(rm));
            
            msubplot(nKeys,1:(coal+1),gs)
            mth = plot(btm,rm,'bO-','linewidth',0.5);
            mstr = 'Motion index (5 min bin)';
        end
    end
    
end

tdstr = '';
tth = nan;
if args.overlay_tdratio
    
    trm = btd_data-nanmin(btd_data);
    trm = datarange*trm/nanmax(trm);
    trm = 100+(trm-nanmean(trm));
    
    msubplot(nKeys,1:(coal+1),gs)
    hold on
    
    tth = plot(bt_td,trm,'rO-','linewidth',1);
    tdstr = 'Th/Del ratio';
end

han = {mth tth};
s = cellfun(@isempty,han);
lstr = {mstr,tdstr};
if any(s)
    leg = legend(han{s},lstr(s),'location','southeast');
    set(leg,'box','off')
end

% Plot theta-delta ratio data
htd = [];

if args.plot_tdratio
    ms = 2;
    i = i+1;
    htd = msubplot(i,1:3,gs);
    plot(t_td,td_data,'O','color',pc,'markersize',ms,'markerfacecolor',pc)
    hold on
    plot([0 0],ylim,'k--')
    grid on
    plot(xlim,[1 1],'r','linewidth',2)
    %     xlim([xl1 xmax])
    ylim([0 quantile(td_data,0.99)])
    ylabel('Theta/Delta ratio')
    box off
end

linkaxes([h htd],'x')

% Plot correlation

if args.plot_corr
    msubplot(nKeys+2,4,gs)
    resp = slopes_pre';
    X2 = btd_data(1:length(resp))';
    plot(X2,resp,'ko','markerfacecolor','none')
    hold on
    lm = fitlm(X2,resp','linear','RobustOpts','on');
    y = predict(lm,X2);
    p = coefTest(lm);
    rs = roundn(lm.Rsquared.Ordinary,-2);
    title(sprintf('%s %s',get_plessthan_str(p),['R^{2} = ' num2str(rs)]),'fontweight','normal','fontsize',10)
    
    plot(X2,y,'r','linewidth',2)
    axis tight
    xlabel('Th/del ratio)')
    ylabel('Norm epsp slope')
    grid on
    box off
    
    
    if args.minfo_exists
        msubplot(nKeys+1,4,gs)
        X1 = (bm_data)';
        lm = fitlm(X1(1:length(resp)),resp,'linear','RobustOpts','on');
        y = predict(lm,X1);
        p = coefTest(lm);
        rs = roundn(lm.Rsquared.Ordinary,-2);
        plot(bm_data(1:length(resp)),resp,'ko','markerfacecolor','none')
        hold on
        plot(X1,y,'r','linewidth',2)
        axis tight
        grid on; box off
        title(sprintf('%s %s',get_plessthan_str(p),['R^{2} = ' num2str(rs)]),'fontweight','normal','fontsize',10)
        
        xlabel('Motion index')
        ylabel('Norm epsp slope')
        
        
        % Plot motion versus th/del ratio
        msubplot(nKeys+2,5,gs)
    
        plot(X1,X2,'ko','markerfacecolor','none')
        hold on
        lm = fitlm(X1,X2,'linear','RobustOpts','on');
        y = predict(lm,X1);
        p = coefTest(lm);
        rs = roundn(lm.Rsquared.Adjusted,-2);
        title(sprintf('%s %s',get_plessthan_str(p),['R^{2} = ' num2str(rs)]),'fontweight','normal','fontsize',10)
        plot(X1,y,'r','linewidth',2)
        axis tight
        xlabel('Motion Index')
        ylabel('Theta/Delta ratio')
    end
end
suptitle(tit_str)
function ph = plot_avg_trace(key,h,col)

axes(h)

[t,y] = fetchn(cstim.FpRespTrace(key),'t','y');
mi = min(cellfun(@length,y));
yi = cellfun(@(x) x(1:mi),y,'uni',false);
yi = [yi{:}]*1000;
t = t{1};
t = t(1:mi)/1000;
ym = median(yi,2);

% Align to trace at 2-3 ms
v = mean(ym(t>=2 & t <=3));
ym = ym-v;
ph = plot(t,ym,'color',col,'linewidth',2);
hold on
xlim([3 30])
yf = ym(t>5 & t<30);
rn = range(yf);
mi = min(yf);
ma = max(yf);
ylf = [mi-rn*0.25 ma+rn*0.25];
xlim([-5 30])
ylim(ylf)
box off
set(gca,'FontSize',12)








