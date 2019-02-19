function plot_nrem_vs_active_epsp(sess_keys,varargin)
% function plot_nrem_vs_active_epsp(sess_keys,varargin)
% MS 2018-05-07

args.motion_th_quantile = 0.5;
args.pre_event_win = 5;
args = parseVarArgs(args,varargin{:});

nSess = length(sess_keys);
leg_str = {'Active','NonREM'};

for iSess = 1:nSess
    % Get the channel keys from this session
    skey = sess_keys(iSess);
    ckeys = fetch(cstim.SlopeEg & cstim.FepspSlope(skey,'smooth_method_num = 0'));
    nKeys = length(ckeys);
    gs = [nKeys,2];
    for iKey = 1:nKeys
        [nra_ratio,nrslopes,actslopes,nrem_trace,act_trace,t] = get_nrem_act_raio(ckeys(iKey),args.pre_event_win,args.motion_th_quantile,1);
        % Bar graph
        msubplot(iKey,1,gs)
        plot_bar_a_vs_b(actslopes,nrslopes,leg_str)
        % Average traces
        msubplot(iKey,2,gs)
        nrh = plot(t,nrem_trace,'r','linewidth',2);
        hold on
        acth = plot(t,act_trace,'k','linewidth',2);
        if iKey ==1
            legend([acth nrh],leg_str,'box','off','location','southeast')
        end
        title(sprintf('nREM/Active ratio: %0.1f',nra_ratio))
        yf = [nrem_trace(:)' act_trace(:)'];
        if ~any(isnan(yf))
        rn = range(yf);
        mi = min(yf);
        ma = max(yf);
        ylf = [mi-rn*0.25 ma+rn*0.25];
        xlim([-5 30])
        ylim(ylf)
        box off
        end
    end
end
ms_suptitle(sprintf('Mouse: %u',skey.animal_id),'yPosition',0.975)