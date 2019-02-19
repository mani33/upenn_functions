function plot_baseline_corr(keys,varargin)
% Plot the correlation between fepsp responses and motion and theta-delta
% ratio for all baselines. Note that no pre-post distinction will be made
% inside the function.
% MS 2017-08-15
args.bin_width = 5;
args.suptitle = '';
args.normalize_motion = 0;
args = parseVarArgs(args,varargin{:});
bwstr = sprintf('slope_bw = %0.2f',args.bin_width);

nKeys = length(keys);
corrdata = struct;

for iKey = 1:nKeys
    key = keys(iKey);
    ckeys = fetch(cstim.SlopeEg & cstim.FepspSlope(key,'smooth_method_num = 0'));
    nChKeys =length(ckeys);
    chdata = struct;
    for iCh = 1:nChKeys
        ckey = ckeys(iCh);
        slopes = fetchn(cstim.SlopeBinned(ckey,bwstr),'y');
        chdata(iCh).sdata = reshape(slopes{:},1,[]);
        x = fetch(cont.MotionBinned(ckey,bwstr),'*');
        if isempty(x)
            x = fetch(cont.MotionInVidBinned(ckey,bwstr),'*');
        end
        chdata(iCh).mdata = x.y(:)';
        tdvals = fetchn(cont.TdrBinned(ckey,bwstr),'y');
        chdata(iCh).tddata = reshape(tdvals{:},1,[]);
    end
    % Get the data from the channel had the best slope-motion correlation.
    bdata = choose_best_correlated(chdata,'motion');
    corrdata.sdata{1,iKey} = normalzeValuesToMean(bdata.sdata);
    if args.normalize_motion
        corrdata.mdata{1,iKey} = normalzeValuesToMean(bdata.mdata);
    else
        corrdata.mdata{1,iKey} = bdata.mdata;
    end
    corrdata.tddata{1,iKey} = bdata.tddata;
end

% Pool across sessions
allslopes = [corrdata.sdata{:}];
allmotion = [corrdata.mdata{:}];
alltddata = [corrdata.tddata{:}];

figure
gs = [3,1];
set(gcf,'Position',[1028,55,483,1293],'color','w')
h = msubplot(1,1,gs);
plot_corr(allmotion,allslopes,h,'Motion Index','Norm epsp slope')

h = msubplot(2,1,gs);
plot_corr(alltddata,allslopes,h,'\theta/\delta ratio','Norm epsp slope')

h = msubplot(3,1,gs);
plot_corr(allmotion,alltddata,h,'Motion Index','\theta/\delta ratio')

ms_suptitle(args.suptitle,'yPosition',0.975)



function plot_corr(X,Y,h,xlab,ylab,titstr)
if nargin < 6
    titstr = '';
end
axes(h)
X = X(:);
Y = Y(:);
lm = fitlm(X,Y,'linear','RobustOpts','on');
y = predict(lm,X);
p = coefTest(lm);
rs = roundn(lm.Rsquared.Ordinary,-2);
plot(X,Y,'k.','markerfacecolor','none')
hold on
plot(X,y,'r','linewidth',2)
axis tight square
grid on; box off
title(sprintf('%s\n%s %s',titstr,get_plessthan_str(p),['R^{2} = ' num2str(rs)]),'fontweight','normal','fontsize',10)
xlabel(xlab)
ylabel(ylab)
xlim(quantile(X(:),[0.01 0.99]))
ylim(quantile(Y(:),[0.01 0.99]))




function nslopes = normalzeValuesToMean(bdata)
nslopes = 100*bdata/mean(bdata);

function [bestChData,bestChInd] = choose_best_correlated(acorr,datatype)

tmp = struct;
nChan = length(acorr);
if nChan == 1
    bestChData = acorr;
    bestChInd = 1;
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
    [~,bestChInd] = max(tmp.rsquared);
    bestChData = acorr(bestChInd);
    disp(['Rsquared values based on ' datatype])
    disp(tmp.rsquared)
    disp('P values')
    disp(tmp.p)
end