function binInd = myhistc(x,edges)

nBins = length(edges)-1;
binInd = zeros(1,length(x));
for i = 1:nBins
    a = edges(i);
    b = edges(i+1);
    bi = x >= a & x < b;
    binInd(bi) = i;
end
