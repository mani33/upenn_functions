function [br,tb] = bin_traces(x,t,bin_width)
% Average the values in x within the given bin_width. Typically used for
% binning epsp responses from field potential recordings.
% Outputs:
% br - binned mean response
% tb - bin center times 

tw = diff(t(1:2));
ns = round(bin_width/tw);
nBins = round(length(x)/ns);
br = cell(1,nBins);
tb = nan(1,nBins);
n = length(x);
for i = 1:nBins
    s1 = (i-1)*ns + 1;
    s2 = min([i*ns n]);
    xs = x(s1:s2);
    xs = cellfun(@(x) x(:),xs,'uni',false);
    br{i} = mean([xs{:}],2);
    tb(i) = mean(t([s1 s2]));
end