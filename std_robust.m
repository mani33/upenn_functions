function s = std_robust(x)
s = median(abs((x-median(x)))/0.6745);