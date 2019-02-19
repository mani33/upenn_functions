function h = plot_ltp_by_ch(keys_pre,keys_post,varargin)

args.plot_his = 1;
args.bin_motion = 1;
args.plot_avg_trace = 1;
args.plot_motion = 1;
args.overlay_motion = 0;
args.plot_tdratio = 1;
args.plot_corr = 1;
args = parseVarArgs(args,varargin{:});

keys_pre = fetch(cstim.SlopeEg & cstim.FepspSlope(keys_pre));
keys_post = fetch(cstim.SlopeEg & cstim.FepspSlope(keys_post));

nKeys = length(keys_pre);
n = length(keys_pre)+double(args.plot_motion)+double(args.plot_tdratio);
binwidth = 5;
havg = nan(1,2);
m = 1;
if args.plot_his
    m = 2;
end

if args.plot_avg_trace
    m = 3;
end
xmax = 180;
ms = 4;
pc = 'k';
plot_raw_resp = 1;
plot_binned = 1;
binwidth = 5;
yl = [25 175];
coal = 2;
gs = [n m+coal];
h = nan(1,nKeys);
for i = 1:nKeys
    %     h = subplot(n,1,i);
    h(i) = msubplot(i,1:(coal+1),gs);
    % Baseline
    key_pre = keys_pre(i);
    key_post = keys_post(i);
    if i == 1
        tit_str = sprintf('M %u: ch %u', key_pre.animal_id,key_pre.chan_num);
    else
        tit_str = sprintf('ch %u', key_pre.chan_num);
    end
    [t,slopes] = fetchn(cstim.FepspSlope(key_pre),'event_ts','fepsp_slope');
    sfac = 100/mean(slopes);
    slopes = slopes * sfac;
    t_minutes = format_time(t);
    ttt = t_minutes-max(t_minutes);
    plot_epsp_resp(ttt,slopes,plot_raw_resp,plot_binned,binwidth,h(i),pc,[],[],[])
    set(gca,'FontSize',12)
    [slopes_pre,~,~] = bin_mouse_data(slopes,ttt,binwidth);
    
    % Post-shock
    [t,slopes] = fetchn(cstim.FepspSlope(key_post),'event_ts','fepsp_slope');
    slopes = slopes * sfac;
    sel = 1:min([180 length(slopes)]);
    if length(sel)~=length(slopes)
        warning('Number of pulses were %u',length(slopes))
    end
    t_minutes = format_time(t);
    t_minutes = t_minutes(sel);
    slopes = slopes(sel);
    plot_epsp_resp(t_minutes,slopes,plot_raw_resp,plot_binned,binwidth,h(i),pc,tit_str,yl,[])
    [slopes_post,~,~] = bin_mouse_data(slopes,t_minutes,binwidth);
    
    plot([0 0],ylim,'k--')
    xl1 = min([-46 ttt(1)]);
    xlim([xl1 xmax])
    grid on
    cylabel('Norm fepsp slope',i==1)
    
    % Plot avg trace
    h2 = msubplot(i,4,gs);
    havg(1) = plot_avg_trace(key_pre,h2,'k');
    havg(2) = plot_avg_trace(key_post,h2,'r');
    if i == 1
        leg = legend(havg,{'Before','After'});
        set(leg,'Box','off','location','southeast')
        xlabel('Time (ms)')
    end
    % Plot histology
    if args.plot_his
        h3 = msubplot(i,5,gs);
        plot_histology(key_pre,1,'r',h3)
    end
    
end

[m_data,tm,td_data,t_td] = get_mouse_motion_and_tdratio_for_ltp(key_pre,key_post);
if args.bin_motion
    ms = 4;
    [bm_data,~,btm] = bin_mouse_data(m_data,tm,binwidth);
    [btd_data,~,bt_td] = bin_mouse_data(td_data,t_td,binwidth);
else
    ms = 2;
end


% Plot motion data
hm = [];
if args.plot_motion
    ms = 2;
    i = i+1;
    chi = i;
    hm = msubplot(i,1:3,gs);
    plot(tm,m_data,'O','color',pc,'markersize',ms,'markerfacecolor',pc)
    hold on
    plot([0 0],ylim,'k--')
    grid on
    linkaxes([h hm],'x')
    xlim([xl1 xmax])
    ylim([0 quantile(m_data,0.99)])
    ylabel('Motion Index')
    box off
end

% Overlay the motion data on the neural responses
if args.bin_motion
    fa = 100;
    trm = btd_data-nanmin(btd_data);
    trm = fa*trm/nanmax(trm);
    trm = 100+(trm-nanmean(trm));
    
    rm = bm_data-nanmin(bm_data);
    rm = fa*rm/nanmax(rm);
    rm = 100+(rm-nanmean(rm));
    
    if args.overlay_motion
    msubplot(nKeys,1:(coal+1),gs)
    mth = plot(btm,rm,'bO-','linewidth',0.5);
    
    
    hold on
    tth = plot(bt_td,trm,'rO-','linewidth',1);
    leg = legend([mth tth],{'Motion index (5 min bin)','Th/Del ratio'},'location','southeast');
    set(leg,'box','off')
    
    %     td_below_one_vals = [];
    %     blah = [];
    % for vikas = 1: length(td_data)
    %     if btd_data(vikas) <= 1
    %         td_below_one_vals(end + 1) = bt_td(vikas);
    %         blah(end + 1) = 50;
    %     else
    %         plot(td_below_one_vals, blah, 'ro', 'LineWidth', 2)
    %         td_below_one_vals = [];
    %         blah = [];
    %     end
    % end
    end
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
    xlim([xl1 xmax])
    ylim([0 quantile(td_data,0.99)])
    ylabel('Theta/Delta ratio')
    box off
end

linkaxes([h htd hm],'x')

% Plot correlation
if args.plot_corr
    msubplot(chi,4,gs)
    resp = [slopes_pre  slopes_post]';
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
    msubplot(chi+1,4,gs)
    
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
    
    % Plot motion versus th/del ratio
    msubplot(chi,5,gs)
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

function ph = plot_avg_trace(key,h,col)

axes(h)

[t,y] = fetchn(cstim.FpRespTrace(key),'t','y');
mi = min(cellfun(@length,y));
yi = cellfun(@(x) x(1:mi),y,'uni',false);
yi = [yi{:}]*1000;
t = t{1};
t = t(1:mi)/1000;
ym = median(yi,2);

% Alighn to trace at 2-3 ms
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



