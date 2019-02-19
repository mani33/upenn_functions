function plot_binned_epsp_resp(ch_key,bin_width,varargin)
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
[t_bin,r_bin,se,t_raw,r_raw] = get_binned_epsp_resp(ch_key,bin_width);
sf = 100/mean(r_raw);
nresp = r_raw * sf;
plot(t_raw([1 end]),[100 100],'r-','linewidth',1)

if args.plot_raw
    plot(t_raw,nresp,'O','color',col,'markersize',ms,'markerfacecolor',col)
    hold on
end

% Get binned responses
r_bin_norm = r_bin * sf;
se_norm = se * sf;
plot(t_bin,r_bin_norm,'O-','color',args.plotCol,'markersize',ms,'markerfacecolor',args.plotCol)
hold on
errorbar(t_bin,r_bin_norm,se_norm,'color',args.plotCol,'linestyle','none')
axis tight

if ~isempty(args.yl)
    ylim(args.yl)
end

title(args.tit_str)
box off

