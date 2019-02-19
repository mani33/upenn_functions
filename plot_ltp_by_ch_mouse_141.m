function h = plot_ltp_by_ch_mouse_141(keys_pre,keys_post1,keys_post2,varargin)

args.plot_his = true;
args.plot_avg_trace = 1;

args = parseVarArgs(args,varargin{:});

keys_pre = fetch(cstim.SlopeEg & cstim.FepspSlope(keys_pre));
keys_post1 = fetch(cstim.SlopeEg & cstim.FepspSlope(keys_post1));
keys_post2 = fetch(cstim.SlopeEg & cstim.FepspSlope(keys_post2));

n = length(keys_pre);
m = 1;
if args.plot_his
    m = 2;
end

if args.plot_avg_trace
    m = 3;
end
pc = 'k';
plot_raw_resp = 1;
plot_binned = 1;
yl = [25 175];

gs = [n m+1];
for i = 1:n
    %     h = subplot(n,1,i);
    h = msubplot(i,1:2,gs);
    % Baseline
    key_pre = keys_pre(i);
    key_post1 = keys_post1(i);
    key_post2 = keys_post2(i);
    
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
    plot_epsp_resp(ttt,slopes,plot_raw_resp,plot_binned,h,pc,[],[],[])
    set(gca,'FontSize',16)
    
    
    % Post-shock
    [t,slopes] = fetchn(cstim.FepspSlope(key_post1),'event_ts','fepsp_slope');
    slopes1 = slopes * sfac;
    sel = 1:100;
    t_minutes = format_time(t);
    t_1 = t_minutes(sel);
    slopes1 = slopes1(sel);
    
    % Get second session for stitching
    [t,slopes2] = fetchn(cstim.FepspSlope(key_post2),'event_ts','fepsp_slope');
    slopes2 = slopes2 * sfac;
    t_2 = format_time(t)+120;
    
    t_minutes = [t_1(:); nan(20,1);t_2(:)];   
    slopes = [slopes1(:); nan(20,1); slopes2(:)];
    plot_epsp_resp(t_minutes,slopes,plot_raw_resp,plot_binned,h,pc,tit_str,yl,[])
    plot([0 0],ylim,'k--')
    xl1 = min([-46 ttt(1)]);
    xlim([xl1 inf])
    % Plot avg trace
    h2 = msubplot(i,3,gs);
    plot_avg_trace(key_pre,h2,'k',[])
    plot_avg_trace(key_post1,h2,'r',key_post2)
    
    % Plot histology
    if args.plot_his
        h3 = msubplot(i,4,gs);
        plot_histology(key_pre,1,'r',h3)
    end
end

function plot_avg_trace(key,h,col,key2)

axes(h)

[t,y] = fetchn(cstim.FpRespTrace(key),'t','y');
if ~isempty(key2)
    [~,y2] = fetchn(cstim.FpRespTrace(key2),'t','y');
    y = [y(1:100); y2];
end
mi = min(cellfun(@length,y));
yi = cellfun(@(x) x(1:mi),y,'uni',false);
yi = [yi{:}]*1000;
t = t{1};
t = t(1:mi)/1000;
ym = mean(yi,2);
% Alighn to trace at 2-3 ms
v = mean(ym(t>=2 & t <=3));
ym = ym-v;
plot(t,ym,'color',col,'linewidth',2)
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
set(gca,'FontSize',16)



