function decFactors = calcDecimationFactors(Fs, cutoffFreq)

targetRate = cutoffFreq * 2 / 0.8;		% 0.8 is the cutoff of the Chebyshev filter in decimate()
coeff = Fs / targetRate;

if (coeff < 2)
    error('Cannot decimate by %g', coeff)
end

% Calculate a series of decimation factors
decFactors = [];
testFactors = 13:-1:2;
while (coeff > 13)
    rems = mod(coeff, testFactors);
    [~, ix] = min(rems);
    decFactors = [decFactors, testFactors(ix)]; %#ok<AGROW>
    coeff = coeff / testFactors(ix);
end

coeff = floor(coeff);
if (coeff >= 2)
    decFactors = [decFactors, coeff];
end


% function out = decimatePackage(data, factor)
% 
% [m,n] = size(data);
% % crop package at a multiple of decimation factor. this is important
% % because otherwise decimate will cause random jitter of up to one sample
% m = fix(m / factor) * factor; 
% out = zeros(m/factor,n);
% for col = 1:n
%     out(:,col) = decimate(data(1:m,col), factor);
% end