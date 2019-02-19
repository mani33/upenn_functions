function [br,se,tb] = bin_resp(x,t,bin_width)
% Average the values in x within the given bin_width. Typically used for
% binning synaptic responses from field potential recordings.
% Outputs:
% br - binned mean response
% se - standard error
% tb - bin center time
tw = diff(t(1:2));
ns = round(bin_width/tw);
nBins = round(length(x)/ns);
br = nan(1,nBins);
se = br;
tb = br;
n = length(x);
for i = 1:nBins
    s1 = (i-1)*ns + 1;
    s2 = min([i*ns n]);
    xs = x(s1:s2);
    br(i) = mean(xs);
    se(i) = std(xs)/sqrt(ns);
    tb(i) = mean(t([s1 s2]));
end