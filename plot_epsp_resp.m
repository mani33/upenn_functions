function plot_epsp_resp(t_minutes,y,plot_raw_resp,plot_binned,bin_width,h,pc,tit_str,yl,fig_size)
% plot_epsp_resp(t_minutes,y,plot_raw_resp,plot_binned,h,pc,tit_str,yl,fig_size)
%

col = [0.8 0.8 0.8];
ms = 6;
subplot(h)
hold on

if plot_raw_resp
    plot(t_minutes,y,'O','color',col,'markersize',ms,'markerfacecolor',col)
    hold on
end
plot(xlim,[100 100],'r-','linewidth',1)
if plot_binned
    [bs,se,tb] = bin_mouse_data(y,t_minutes,bin_width);
    plot(tb,bs,'O','color',pc,'markersize',ms,'markerfacecolor',pc)
    hold on
    errorbar(tb,bs,se,'color',pc,'linestyle','none')
end
axis tight
% plot(xlim,[100 100],'k--')
if ~isempty(yl)
    ylim(yl)
end

title(tit_str)
box off
if ~isempty(fig_size)
    set_fig_size(gcf,fig_size)
end