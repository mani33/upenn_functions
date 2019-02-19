function plot_ltp_data_by_ch(key_pre,key_post,key_nextday_baseline,varargin)
% function plot_ltp_data_by_ch(key_pre,key_post,key_nextday_baseline,varargin)
% 2017-08-11, MS
% 2018-04-25, MS

args.motion_idx_quantiles = [0 0.95];
args.motion_th_quantile = 0.5;
args.use_all_data = false;
args.max_ch = inf;
args.invert_time = 0;
args.rec_gap = 5; % gap between pre and post session
args.overlay_motion = 0;
args.overlay_tdr = 0;
args.yl = [50 175];
args.plot_corr = 1;
args.plot_his = 0;
args.tit = '';
args.plot_errorbar = 0;
args.plot_best_ch = 0;
args.plot_all_chan_tdr = 0;
args.show_pval = 0;
args.save = false;
args.save_dir = '';
args.ch_sel_method = 'ltp_magnitude';
args.keep_just_motion = false;
args.slope_bw = 5;
args.use_pre_event_motion = 1;
args.use_pre_event_tdr = 1;
args.pre_event_win = 5;
args = parseVarArgs(args,varargin{:});

ckeys_pre = fetch(cstim.SlopeEg & cstim.FepspSlope(key_pre,'smooth_method_num = 0'));
if ~isempty(key_post)
    ckeys_post = fetch(cstim.SlopeEg & cstim.FepspSlope(key_post,'smooth_method_num = 0'));
    if ~isempty(key_nextday_baseline)
        ckeys_day2 = fetch(cstim.SlopeEg & cstim.FepspSlope(key_nextday_baseline,'smooth_method_num = 0'));
    else
        ckeys_day2 = [];
    end
    invtime = true;
    if args.plot_best_ch
        [~,ckeys_pre,ckeys_post,ckeys_day2] = get_best_channel(key_pre,key_post,ckeys_day2,args.ch_sel_method,args);
        args.max_ch = 1;
    end
else
    ckeys_post = [];
    ckeys_day2 = [];
    invtime = false;
end

figure
set(gcf,'Position',[1991,311,1592,956])
nChKeys = min(length(ckeys_pre),args.max_ch);

% Plotting related
respCol = 5;
traceNhisCol = 2 + double(args.plot_his);
tdrCol = 0;
if args.plot_all_chan_tdr
    tdrCol = 6;
end
npcol = respCol + traceNhisCol + tdrCol;
nprow = nChKeys + 4;
day2plot = double(~isempty(key_nextday_baseline))*2;
gs = [nprow npcol+day2plot];

% raw slope color/marker size
sraw.col = [0.9 0.9 1];
sraw.ms = 3;
sbin.col = 'k';
sbin.ms = 3;
precol = 'k';
postcol = 'r';
day2col = [0.2 0.2 1];
% Motion
mraw.col = 'k';
mraw.ms = 1;
mbin.col = 'k';
mbin.ms = 1;
if nChKeys < 5
    mraw.ms = 2;
    mbin.ms = 2;
end
% TD ratio
tdraw.col = 'k';
tdraw.ms = 1;
tdbin.col = 'k';
tdbin.ms = 1;
if nChKeys < 5
    tdraw.ms = 2;
    tdbin.ms = 2;
end
% Overlay
olaym.col = [0 0.5 0];
olaym.ms = 2;
olaytd.col = 'b';
olaytd.ms = 2;
olaySpread = 70;

% Correlation related
adata = struct;
% For choosing the best tdratio channel
atddata = struct;

hh = nan(1,1);
if ~isempty(ckeys_post)
args.motion_idx_bound = get_common_motion_idx_bound(ckeys_pre(1),ckeys_post(1),args);
else
    args.motion_idx_bound = get_common_motion_idx_bound(ckeys_pre(1),[],args);
end

for iCh = 1:nChKeys
    ckey_pre = ckeys_pre(iCh);
    tit_str = sprintf('ch %u', ckey_pre.chan_num);

    h = nan;
    avgStr = cell(1,1);
    
    [sdata1, mdata1, tddata1,pre_ev_mov_idx1,pre_ev_tdr1] = get_all_data(ckey_pre,args,'invert_time',invtime);
    sfr = 100/nanmean(sdata1.raw_val);
    sfb = 100/nanmean(sdata1.bin_val);
    
    % Plot fepsp responses
    hh(iCh) = msubplot(iCh,1:respCol,gs);
    plot_slopes(sdata1,sraw,sbin,sfr,sfb,args.plot_errorbar)
    xmin = sdata1.raw_t(1);
    xmax = sdata1.raw_t(end);
    
    % Plot avarage response
    xh = msubplot(iCh,respCol+day2plot+(1:2),gs);
    h(1) = plot_avg_trace(ckey_pre,xh,precol,args,pre_ev_mov_idx1);
    cxlabel('Time (ms)',iCh==1)
    avgStr{1} = 'pre';
    % Post session
    sdata2.bin_val = [];
    sdata2.bin_t = [];
    
    mdata2.bin_val = [];
    mdata2.bin_t = [];
    tddata2 = struct('bin_val',[],'bin_t',[],'se',[],'raw_val',[],'raw_t',[],'th',[],'del',[],...
        'se_th',[],'se_del',[]);
    if ~isempty(ckeys_post)
        hh(iCh) = msubplot(iCh,1:respCol,gs);
        ckey_post = ckeys_post(iCh);
        [sdata2, mdata2, tddata2, pre_ev_mov_idx2,pre_ev_tdr2] = get_all_data(ckey_post,args);
        plot_slopes(sdata2,sraw,sbin,sfr,sfb,args.plot_errorbar)
        xmax = sdata2.raw_t(end);
        % Plot average response
        h(2) = plot_avg_trace(ckey_post,xh,postcol,args,pre_ev_mov_idx2);
        avgStr{2} = 'post';
    end
    
    %% Next day baseline
    if ~isempty(ckeys_day2)
        % Post session
        sdata3.bin_val = [];
        sdata3.bin_t = [];
        
        mdata3.bin_val = [];
        mdata3.bin_t = [];
        tddata3 = struct('bin_val',[],'bin_t',[],'se',[],'raw_val',[],'raw_t',[],'th',[],'del',[],...
            'se_th',[],'se_del',[]);
        if ~isempty(ckeys_day2)
            hh(iCh) = msubplot(iCh,(1:2)+respCol,gs);
            ckey_d2 = ckeys_day2(iCh);
            [sdata3, mdata3, tddata3,pre_ev_mov_idx3,pre_ev_tdr3] = get_all_data(ckey_d2,args);
            plot_slopes(sdata3,sraw,sbin,sfr,sfb,args.plot_errorbar)
            xlim([0 sdata3.bin_t(end)+diff(sdata3.bin_t(1:2))])
            ylim(args.yl)
            box off
            plot(xlim,[100 100],'r','linewidth',1)
            % Plot average response
            h(3) = plot_avg_trace(ckey_d2,xh,day2col,args,pre_ev_mov_idx2);
            avgStr{3} = 'day2';
        end
    end
    
    if iCh==1, leg = legend(h,avgStr); set(leg, 'box', 'off','location','southeast');end
    hh(iCh) = msubplot(iCh,1:respCol,gs);
    xlim(hh(iCh),[xmin xmax])
    ylim(args.yl)
    plot(xlim,[100 100],'r','linewidth',1)
    title(tit_str)
    cylabel('Normalized fepsp slope',iCh==1)
    box off
    grid off
    plot([0 0],ylim,'r--')
    xli = xlim;
    % Save data for correlations
    adata(iCh).sdata = sfb*[sdata1.bin_val sdata2.bin_val]';
    adata(iCh).mdata = [mdata1.bin_val mdata2.bin_val]';
    adata(iCh).tddata = [tddata1.bin_val tddata2.bin_val]';
    adata(iCh).t = [sdata1.bin_t sdata2.bin_t]';
    
    % Save tdratio data for plotting the 'best' one later
    atddata(iCh).tddata1 = tddata1;
    atddata(iCh).tddata2 = tddata2;
    
    % Plot histology
    if args.plot_his
        % Recording electrode
        hStim = msubplot(iCh,respCol+4,gs);
        plot_histology(ckey_pre,1,'r',hStim,'Rec')
        % Stimulating electrode
        hStim = msubplot(1,respCol+2,gs);
        plot_histology(ckey_pre,0,'r',hStim,'Stim','stim',1)
    end
end

% Choose the best channel that shows the most correlation with either
% motion or tdratio.
% best_chan_ind - bci
if args.plot_corr
    [corr,bci_m] = choose_best_correlated(adata,'motion');
    [~,bci_td] = choose_best_correlated(adata,'tdratio');
    
    tddata1 = atddata(bci_td).tddata1;
    tddata2 = atddata(bci_td).tddata2;
end
%% Motion data

if ~isempty(mdata1.bin_val)
    iCh = iCh + 1;
    hh(iCh) = msubplot(iCh,1:respCol,gs);
    % Pre
    plot(mdata1.raw_t, mdata1.raw_val,'O','color',mraw.col,'markersize',mraw.ms,'markerfacecolor',mraw.col)
    hold on
   
    % Post
    if ~isempty(ckeys_post)&& ~isempty(mdata2)
        plot(mdata2.raw_t, mdata2.raw_val,'O','color',mraw.col,'markersize',mraw.ms,'markerfacecolor',mraw.col)
    end
    mtop = quantile(mdata1.raw_val,0.90);
    ylim([-0.05*mtop mtop])
    ylabel('Motion Index','color',olaym.col)
    box off
    xlim(xli)
    
    % Motion thresholded data
    iCh = iCh + 1;
    msubplot(iCh,1:respCol,gs);
    mv = pre_ev_mov_idx1;
    mv(mv < args.motion_idx_bound(1) | mv > args.motion_idx_bound(2)) = nan;
    plot(sdata1.raw_t,pre_ev_mov_idx1,'k.','markersize',12)
    hold on
    plot(sdata1.raw_t,mv,'r.','markersize',12)
    ylabel('Pre-stim motion idx')
    if ~isempty(ckeys_post)&& ~isempty(mdata2)
        % Motion thresholded data
        mv = pre_ev_mov_idx2;
        mv(mv < args.motion_idx_bound(1) | mv > args.motion_idx_bound(2)) = nan;
        plot(sdata2.raw_t,pre_ev_mov_idx2,'k.','markersize',12)
        hold on
        plot(sdata2.raw_t,mv,'r.','markersize',12)
        mtop = quantile(pre_ev_mov_idx2,0.95);
        ylim([-0.05*mtop mtop])
        box off
        xlim(xli)
    end
    
    % Day2
    if ~isempty(ckeys_day2) && ~isempty(mdata3)
        msubplot(iCh,(1:2)+respCol,gs);
        plot(mdata3.raw_t, mdata3.raw_val,'O','color',mraw.col,'markersize',mraw.ms,'markerfacecolor',mraw.col)
        box off
        axis tight
        ylim([-0.05*mtop mtop])
     
        % Motion thresholded data
        mv = pre_ev_mov_idx3;
        mv(mv < args.motion_idx_bound(1) | mv > args.motion_idx_bound(2)) = nan;
        plot(sdata3.raw_t,pre_ev_mov_idx3,'k.','markersize',12)
        hold on
        plot(sdata3.raw_t,mv,'r.','markersize',12)
        mtop = quantile(pre_ev_mov_idx3,0.95);
        ylim([-0.05*mtop mtop])
        box off
        xlim(xli)
    end
    % Overlay motion
    if args.overlay_motion
        msubplot(bci_m,1:respCol,gs)
        plot(mdata1.bin_t, scale_shift(mdata1.bin_val,olaySpread), 'O-','color',olaym.col,'markersize',olaym.ms,'markerfacecolor',olaym.col)
        if ~isempty(ckeys_post)&& ~isempty(mdata2)
            plot(mdata2.bin_t, scale_shift(mdata2.bin_val,olaySpread), 'O-','color',olaym.col,'markersize',olaym.ms,'markerfacecolor',olaym.col)
        end
    end
end
%% Plot td ratio for the last chan
iCh = iCh + 1;
hh(iCh) = msubplot(iCh,1:respCol,gs);
% Pre
plot(tddata1.raw_t, tddata1.raw_val,'O','color',tdraw.col,'markersize',tdraw.ms,'markerfacecolor',tdraw.col)
hold on

% Post

if ~isempty(ckeys_post)
    plot(tddata2.raw_t, tddata2.raw_val,'O','color',tdraw.col,'markersize',tdraw.ms,'markerfacecolor',tdraw.col)
end
ylim([0 quantile(tddata1.raw_val,0.99)])
plot(xlim,[1 1],'r','linewidth',1)
ylabel('\theta/\delta ratio','color',olaytd.col)
xlabel('Time (min)')
box off
% linkaxes(hh,'x')
xlim(xli)

set(gcf,'color','w')

% Day2
if ~isempty(ckeys_day2) && ~isempty(tddata3)
    msubplot(iCh,(1:2)+respCol,gs);
    plot(tddata3.raw_t, tddata3.raw_val,'O','color',tdraw.col,'markersize',tdraw.ms,'markerfacecolor',tdraw.col)
    %         ylabel('Motion Index','color',olaym.col)
    hold on
    plot(xlim,[1 1],'r','linewidth',1)
    
    box off
    xlim([0 tddata3.bin_t(end)+diff(tddata3.bin_t(1:2))])
    ylim([0 quantile(tddata1.raw_val,0.99)])
end

%% Overlay of binned td ratio
if args.overlay_tdr
    msubplot(bci_td,1:respCol,gs)
    v1 = scale_shift(tddata1.bin_val,olaySpread);
    plot(tddata1.bin_t,v1 , 'O-','color',olaytd.col,'markersize',olaytd.ms,'markerfacecolor',olaytd.col)
    if ~isempty(ckeys_post)&& ~isempty(mdata2)
        v2 = scale_shift(tddata2.bin_val,olaySpread);
        plot(tddata2.bin_t, v2, 'O-','color',olaytd.col,'markersize',olaytd.ms,'markerfacecolor',olaytd.col)
    end
end

% Plot pre-event tdr
iCh = iCh + 1;
msubplot(iCh,1:respCol,gs);
tv = pre_ev_tdr1;
mv = pre_ev_mov_idx1;
tv(mv < args.motion_idx_bound(1) | mv > args.motion_idx_bound(2)) = nan;
plot(sdata1.raw_t,pre_ev_tdr1,'k.','markersize',12)
hold on
plot(sdata1.raw_t,tv,'r.','markersize',12)
if ~isempty(ckeys_post)&& ~isempty(mdata2)
    % Motion thresholded data
    mv = pre_ev_mov_idx2;
    tv = pre_ev_tdr2;
    tv(mv < args.motion_idx_bound(1) | mv > args.motion_idx_bound(2)) = nan;
    plot(sdata2.raw_t,pre_ev_tdr2,'k.','markersize',12)
    hold on
    plot(sdata2.raw_t,tv,'r.','markersize',12)
    mtop = quantile([pre_ev_tdr1(:)' pre_ev_tdr2(:)'],0.95);
    ylim([-0.05*mtop mtop])
    box off
    xlim(xli)
   ylabel('Pre stim \theta/\delta ratio','color',olaytd.col)
end


%%
if args.plot_corr
    %     msubplot(nChKeys+1,5,gs)
    
    X1 = corr.tddata;
    Y = corr.sdata;
    %     % TD ratio versus slopes
    %     plot(X1,Y,'k.','markerfacecolor','none')
    %     hold on
    %     lm = fitlm(X1,Y,'linear','RobustOpts','on');
    %     y = predict(lm,X1);
    %     p = coefTest(lm);
    %     rs = roundn(lm.Rsquared.Ordinary,-2);
    %     title(sprintf('%s %s',get_plessthan_str(p),['R^{2} = ' num2str(rs)]),'fontweight','normal','fontsize',10)
    %
    %     plot(X1,y,'r','linewidth',2)
    %     axis tight square
    %     xlabel('\theta/\delta ratio')
    %     ylabel('Norm epsp slope')
    %     grid on
    %     box off
    
    if ~isempty(corr.mdata)
        X2 = corr.mdata;
        msubplot(nChKeys+1,respCol+day2plot+(1:2),gs);
        lm = fitlm(X2,Y,'linear','RobustOpts','on');
        y = predict(lm,X2);
        p = coefTest(lm);
        rs = roundn(lm.Rsquared.Ordinary,-2);
        plot(X2,Y,'k.','markerfacecolor','none')
        hold on
        plot(X2,y,'r','linewidth',2)
        axis tight square
        grid on; box off
        if args.show_pval
            title(sprintf('%s %s',get_plessthan_str(p),['R^{2} = ' num2str(rs)]),'fontweight','normal','fontsize',10)
        end
        xlabel('Motion index')
        ylabel('Norm epsp slope')
        
        %         % Plot motion versus th/del ratio
        %         msubplot(nChKeys+1,6,gs)
        %
        %         plot(X2,X1,'k.','markerfacecolor','none')
        %         hold on
        %         lm = fitlm(X2,X1,'linear','RobustOpts','on');
        %         y = predict(lm,X2);
        %         p = coefTest(lm);
        %         rs = roundn(lm.Rsquared.Adjusted,-2);
        %         title(sprintf('%s %s',get_plessthan_str(p),['R^{2} = ' num2str(rs)]),'fontweight','normal','fontsize',10)
        %         plot(X2,y,'r','linewidth',2)
        %         axis tight square
        %         xlabel('Motion Index')
        %         ylabel('\theta/\delta ratio')
        %         box off
    end
end


%% Plot individual channel theta-delta data
if args.plot_all_chan_tdr
    tmp = struct;
    for iCh = 1:nChKeys
        tdpre = atddata(iCh).tddata1;
        tdpost = atddata(iCh).tddata2;
        t = [tdpre.bin_t nan tdpost.bin_t];
        th = [tdpre.th nan tdpost.th];
        se_th = [tdpre.se_th nan tdpost.se_th];
        
        del = [tdpre.del nan tdpost.del];
        se_del = [tdpre.se_del nan tdpost.se_del];
        
        tdr = [tdpre.bin_val nan tdpost.bin_val];
        se_tdr = [tdpre.se nan tdpost.se];
        
        % Single channel plot
        % Theta
        tmp.hth(iCh) = msubplot(iCh,7:8, gs);
        plot(t,th,'k.--');
        hold on
        %     axis tight
        %     errorbar(t,th,se_th,'k','linestyle','none')
        ctitle('\theta (4-8Hz) band power',iCh==1)
        %     ylim([0 8]*1e-8)
        % Delta
        tmp.hdel(iCh) = msubplot(iCh,9:10, gs);
        plot(t,del,'k.--');
        hold on
        %     errorbar(t,del,se_del,'k','linestyle','none')
        ctitle('\delta (1-4Hz) band power',iCh==1)
        box off
        %     axis tight
        %     ylim([0 8]*1e-7)
        % Ratio
        tmp.htdr(iCh) = msubplot(iCh,11:12,gs) ;
        plot(t,tdr,'k.--');
        hold on
        %     errorbar(t,tdr,se_tdr,'color','k','linestyle','none')
        ctitle('\theta/\delta ratio',iCh==1)
        %     leg = legend([hth hdel htdr],{'\theta','\delta','\theta/\delta'});
        %     set(leg,'box','off')
        %     ylim([0 1.25])
        box off
        
        linkaxes([tmp.hth(iCh) tmp.hdel(iCh) tmp.htdr(iCh)],'x')
        xlim([t(1) max(t)])
    end
end

if isempty(ckeys_post)
    post_str = [];
else
    [~,pstr] = fileparts(fetch1(acq.Sessions(key_post),'session_path'));
    pstr = strrep(pstr,'_','\_');
    post_str = [pstr ' (post)'];
end

[~,pre_str] = fileparts(fetch1(acq.Sessions(key_pre),'session_path'));
pre_str = strrep(pre_str,'_','\_');
% stit = [sprintf('Mouse # %u:      ',key_pre.animal_id) [pre_str ' (pre)  ' post_str] '  ' args.tit];
if isempty(ckeys_day2)
    d2str = [];
else
    [~,day2_str] = fileparts(fetch1(acq.Sessions(key_nextday_baseline),'session_path'));
    day2_str = strrep(day2_str,'_','\_');
    d2str = [day2_str ' (day2) '];
end
stit = [sprintf('Mouse # %u:      ',key_pre.animal_id) [pre_str ' (pre)  ' post_str '  ' d2str] '  ' args.tit];

suptitle(stit)
if args.save
    fn = sprintf('mouse_%u_%s',key_pre.animal_id,args.tit);
    saveas(gcf,fn,'jpeg')
end

function [corr,best] = choose_best_correlated(acorr,datatype)

tmp = struct;
nChan = length(acorr);
if nChan == 1
    corr = acorr;
    best = 1;
else
    for iChan = 1:nChan
        dd = acorr(iChan);
        Y = dd.sdata;
        switch datatype
            case 'motion'
                X = dd.mdata;
            case 'tdratio'
                X = dd.tddata;
        end
        lm = fitlm(X,Y,'linear','RobustOpts','on');
        tmp.rsquared(iChan) = lm.Rsquared.Ordinary;
        tmp.p(iChan) = coefTest(lm);
    end
    [~,best] = max(tmp.rsquared);
    corr = acorr(best);
    disp(['Rsquared values based on ' datatype])
    disp(tmp.rsquared)
    disp('P values')
    disp(tmp.p)
end


function v = scale_shift(v,olaySpread)
% Scale and shift
v = v - nanmin(v);
v = olaySpread*v/nanmax(v);
v = 100+(v-nanmean(v));

function plot_slopes(sdata,sraw,sbin,sfr,sfb,plot_eb)

% Raw responses
rs = sdata.raw_val*sfr;
plot(sdata.raw_t, rs,'O','color',sraw.col,'markersize',sraw.ms,'markerfacecolor',sraw.col)
hold on
% Binned
brs = sdata.bin_val*sfb;
plot(sdata.bin_t,brs ,'O','color',sbin.col,'markersize',sbin.ms,'markerfacecolor',sbin.col,'linewidth',1.5)
if plot_eb
    % errorbar
    errorbar(sdata.bin_t, brs, sdata.se*sfb,'color',sbin.col,'linestyle','none')
end

%
function [sdata, mdata, tddata, pre_ev_mov_idx, pre_ev_tdr] = get_all_data(ch_key,args,varargin)

% Get raw and binned slopes
[sdata,pre_ev_mov_idx] = get_binned_epsp_resp(ch_key,args,varargin{:});
% Get motion data if exists
mdata = get_binned_motion(ch_key,args,varargin{:});
% Get theta delta ratio
[tddata ,pre_ev_tdr] = get_binned_tdrdata(ch_key,args,varargin{:});


function ph = plot_avg_trace(key,h,col,args,pre_ev_mov_idx,varargin)

args.linewidth = 1.5;
args.linestyle = '-';
args = parseVarArgs(args,varargin{:});

axes(h)

[t,y,et] = fetchn(cstim.FpRespTrace(key),'t','y','event_ts');
if key.animal_id == 141 && double(key.session_start_time) == 5333866070
    nBad = 20;
    filler_raw = cell(nBad,1);
    t((end-(nBad-1)):end) = [];
    y((end-(nBad-1)):end) = [];
    et((end-(nBad-1)):end) = [];
    % second session to stitch
    key2 = fetch(acq.Ephys(get_key_from_sess_ts('2017-07-17_15-39-25')));
    mch_key = key2;
    mch_key.chan_num = key.chan_num;
    mch_key.slope_win = key.slope_win;
    mch_key.smooth_method_num = key.smooth_method_num;
    [t2,y2,et2] = fetchn(cstim.FpRespTrace(mch_key),'t','y','event_ts');
    t = [t; filler_raw;t2];
    y = [y; filler_raw;y2];
    et = [et; nan(nBad,1); et2];
end
sel = pre_ev_mov_idx > args.motion_idx_bound(1) & pre_ev_mov_idx < args.motion_idx_bound(2);
y = y(sel);
et = et(sel);
t = t(sel);
t = t{1};
et = double(et);
if args.keep_just_motion
    % Apply motion based selection of slopes
    d = fetch(cont.MovSegManual(key),'*');
    if isempty(d)
        d = fetch(cont.MovSegInVidManual(key),'*');
    end
    nSeg = length(d.seg_begin);
    % Start with a NaN matrix
    tmp.val = cell(size(et));
    if ~isnan(d.seg_begin(1))
        for iSeg = 1:nSeg
            sel = et >= d.seg_begin(iSeg) & et <= d.seg_end(iSeg);
            tmp.val(sel) = y(sel);
        end
    end
    y = tmp.val;
end

ey = cellfun(@isempty,y);
y = y(~ey);
mi = min(cellfun(@length,y));
yi = cellfun(@(x) x(1:mi),y,'uni',false);
yi = [yi{:}]*1000;

t = t(1:mi)/1000;
ym = median(yi,2);

% Align to trace at 2-3 ms
v = mean(ym(t>=2 & t <=3));
ym = ym-v;
ph = plot(t,ym,'color',col,'linewidth',args.linewidth,'linestyle',args.linestyle);
hold on
xlim([3 30])
yf = ym(t>5 & t<30);
rn = range(yf);
mi = min(yf);
ma = max(yf);
ylf = [mi-rn*0.25 ma+rn*0.25];
xlim([-5 30])
if key.animal_id == 202
    ylim([-0.1236    1.3361])
else
    ylim(ylf)
end
box off
set(gca,'FontSize',8)
ylabel('mV')