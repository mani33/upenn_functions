function varargout = plot_epsp_vs_motion_tdr(ch_key,xaxis_var,yaxis_var,varargin)
% function h = plot_epsp_vs_motion_tdr(ch_key,xaxis_var,yaxis_var)
% Plot motion or theta/delta ratio against epsp slope or popspike for a
% given channel key
% MS 2018-05-02
%
args.rel_history_begin = -5; 
args.rel_history_end = 0; 
args.motion_th_quantile = 0;
args.axes = [];
args.xscale = 'linear';
args.col = 'k';
args.fit = true;
args = parseVarArgs(args,varargin{:});
if isempty(args.axes)
    args.axes = gca;
end

switch yaxis_var
    case 'slope'
        [y,yt] = fetchn(cstim.FepspSlope(ch_key),'fepsp_slope','event_ts');
        y = -y;
    case 'popspike'
        [y,yt] = fetchn(cstim.FepspSlope(ch_key),'popspike_height','event_ts');
    otherwise
        error('Only slope or popspike allowed')
end
yt = double(yt);
switch xaxis_var
    case 'motion'
        xd = fetch(cont.PreEventMotion(ch_key),'dist_var','event_ts');
        if isempty(xd)
            xd = fetch(cont.MotionInVid(ch_key),'dist_var','t');
        end
        if isempty(xd)
            error('Motion or MotionInVid table tuple does not exist')
        end
        xv = [xd.dist_var];
        q = quantile(xv,args.motion_th_quantile);
        xv(xv < q) = nan;
        xt = double([xd.event_ts]);
    case 'tdratio'
        xd = fetch(cont.TDratio(ch_key),'td_ratio','t');
        if isempty(xd)
            error('TDratio table tuple does not exist')
        end
        xv = xd.td_ratio;
        xt = xd.t;
    otherwise
        error('only motion or tdratio is allowed')
end

x = get_x_vals(xv,xt,yt,args);
axes(args.axes)
if strcmp(args.xscale,'log')
    x = log10(x);
end
h = plot(x,y,'O','color',args.col);
hold on

if args.fit
    x(isinf(x)) = nan;
    lm = fitlm(x,y,'linear','RobustOpts','on');
    yhat = predict(lm,x);
    p = coefTest(lm);
    rs = roundn(lm.Rsquared.Ordinary,-2);
    plot(x,yhat,'r','linewidth',2)
    legend(sprintf('%0.2f',rs))
    disp(p)
    xlabel('Motion')
    ylabel('Epsp response (slope)')
end

if nargout
    varargout{1} = h;
end

function x = get_x_vals(xv,xt,yt,args)
% Here we compute the average value of x (motion or tdratio) in a short
% window that immediately precedes the stim pulse(yt).

nPulse = length(yt);
x = nan(nPulse,1);
for iPulse = 1:nPulse
    ct = yt(iPulse);
    twin = [args.rel_history_begin args.rel_history_end]*1e6 + ct;
    x(iPulse) = mean(xv(xt >= twin(1) & xt <= twin(2)));
end
if any(isnan(x))
    fprintf('Number of pulses without motion/tdr data = %u out of %u\n',length(find(isnan(x))),nPulse)
end



