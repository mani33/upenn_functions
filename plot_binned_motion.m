function plot_binned_motion(sesskey,bin_width,varargin)
% plot_binned_epsp_resp(ch_key,bin_width,varargin)
%
args.plot_raw = 1;
args.h = [];
args.pc = 'k';
args.tit_str = '';
args.yl = []; % ylimit
args.plotCol = 'k';
args.invert_time = 0;
args = parseVarArgs(args,varargin{:});

tfac = (1e-6)/60; % to convert to min
col = [0.8 0.8 0.8];
ms = 6;
if isempty(args.h)
    h = gca;
end
subplot(h)
hold on
rd = fetch(cstim.FepspSlope(ch_key,'smooth_method_num = 0'),'event_ts','fepsp_slope');
traw = double([rd.event_ts]);
t0 = traw(1);
if args.invert_time
    t0 = traw(end);
end
rresp = [rd.fepsp_slope];
sf = 100/mean(rresp);
nresp = rresp * sf;
t_nobin = (traw-t0)*tfac; % in minutes

plot(t_nobin([1 end]),[100 100],'r-','linewidth',1)

if args.plot_raw
    plot(t_nobin,nresp,'O','color',col,'markersize',ms,'markerfacecolor',col)
    hold on
end

% Get binned responses
bwstr = sprintf('slope_bw = %0.2f',bin_width);
rd = fetch(cstim.SlopeBinned(ch_key,'smooth_method_num = 0',bwstr),'t','y','se');
if isempty(rd)
    warning('no SlopeBinned data found for given key')
    return
end
rd.t = double(rd.t);
t_bin = (rd.t-t0)*tfac; % in minutes
plot(t_bin,rd.y*sf,'O','color',args.plotCol,'markersize',ms,'markerfacecolor',args.plotCol)
hold on
errorbar(t_bin,rd.y*sf,rd.se*sf,'color',args.plotCol,'linestyle','none')
axis tight

if ~isempty(args.yl)
    ylim(args.yl)
end

title(args.tit_str)
box off

