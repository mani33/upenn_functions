function p = test_baseline_stability_kolmogorov_dimitri(d,varargin)
% function p = test_baseline_stability_kolmogorov_dimitri(d)
% Test baseline stability of field recording normalized slope values based
% on the idea given by Dimitri Yatsenko.
% 2018-02-27

args.debug = true;
args.keep_just_motion = false;
args = parseVarArgs(args,varargin{:});
%% First remove NaN's
nanInd = (isnan(d));
d = d(~nanInd);
nNan = length(find(nanInd));
if nNan>0
    warning('%0.0d NaN values were removed',nNan)
end
n = length(d);
nBoot = 10000;
model_mean = 100;
assert(abs(model_mean-mean(d)) < 10, 'You need to normalize the slopes to 100%')
csdata = cumsum(d);
% csdata = csdata/csdata(end);
cso = csdata;
csmodel = cumsum(model_mean*ones(size(d)));
% csmodel = csmodel/csmodel(end);
D = max(abs(csmodel-csdata));
Di = nan(1,nBoot);
for i = 1:nBoot
    sel = ceil(rand(1,n)*n);
    ds = d(sel);
%     if i==1
%     figure(1222)
%     plot(ds,'kO')
%     hold on
%     plot(d,'r*')
%     plot(1:length(d),100*ones(size(d)),'b.','markersize',5)
%     hold off
%     end
    csdata = cumsum(ds);
    Di(i) = max(abs(csdata-csmodel));
end
if args.debug
    figure
    subplot(1,2,1)
    plot(cso,'r.-','markersize',5);
    hold on
    plot(csmodel,'.','markerfacecolor','b','markersize',5)
    xlabel('Time (min)'); ylabel('Cumulative sum of slopes')
    subplot(1,2,2)
    hist(Di,50)
    hold on
    plot([D D],ylim,'r')
    xlabel('D statistic'); ylabel('Freq')
end
p = length(find(Di>D))/nBoot;
