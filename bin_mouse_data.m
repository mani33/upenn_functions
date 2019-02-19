function [bs,se,tb] = bin_mouse_data(x,t,b)
% [binnedSlope, std_error, time] = bin_mouse_data(slopes,t_minutes,bin_width_min)

tw = diff(t(1:2));
ns = round(b/tw);
sz = length(x);
nBins = round(sz/ns);
bs = nan(1,nBins);
se = bs;
tb = bs;

for i = 1:nBins
    s1 = (i-1)*ns + 1;
    s2 = i*ns;
    if s2 > sz
        s2 = sz;
    end
    xs = x(s1:s2);
    bs(i) = mean(xs);
    se(i) = std(xs)/sqrt(ns);
    tb(i) = t(s1);
end
tb = tb+b/2;





