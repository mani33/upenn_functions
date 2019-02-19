function plot_bar_a_vs_b(adata,bdata,condStr,varargin)
args.col_a = 'k';
args.col_b = 'k';
args.ylabel = '';
args.show_stat = false;
args.title = '';
args.boxplot = 0;
args.paired = false;
args.show_legend = false;
args.axlim = [0 3 -inf inf];
args.face_col1 = [0.8 0.8 0.8];
args.face_col2 = [1 0 0];

args = parseVarArgs(args,varargin{:});

lena = length(adata);
lenb = length(bdata);
m = [nanmean(adata) nanmean(bdata)];
se = [nanstd(adata) nanstd(bdata)]./[sqrt(lena) sqrt(lenb)];
if args.boxplot
    boxplot([adata bdata],'labels',condStr,'notch','off')
else
    x = bar(1,m(1),args.col_a);
    set(x,'facecolor',args.face_col1,'linewidth',1)
    hold on
    
    y = bar(2,m(2),args.col_b);
    set(y,'facecolor',args.face_col2,'linewidth',1);
    errorbar([1,2],m,se,'linestyle','none','color','k')
    
    axis(args.axlim)
    
    if args.paired
        X = repmat([1 2],length(adata),1);
        Y = [adata bdata];
        plot(X',Y','kO-','markerfacecolor','k')
    else
        plot(ones(1,lena),adata,'O','color','k','markerfacecolor','k')
        plot(2*ones(1,lenb),bdata,'kO','markerfacecolor','k')
    end
end
hold on
if args.paired
        X = repmat([1 2],length(adata),1);
        Y = [adata bdata];
        plot(X',Y','kO-','markerfacecolor','k')
end
set(gca,'linewidth',1,'fontsize',14)
box off
if args.show_legend
    leg = legend([x,y],condStr);
    set(leg,'box','off')
end
if args.show_stat
    [h,p] = ttest2(adata,bdata);
    text(2.2,88,sprintf('n = %0.0f,%0.0f;  p = %0.2f',lena,lenb,p))
end
ylabel(args.ylabel)
set(gca,'xtick',[1 2],'xticklabel',condStr)
title(args.title)
