function plot_freezing_vs_ltp(keys,ch)

tmp = struct;
n = length(keys);
for i = 1:n
   key = keys(i);
   slopes = fetchn(cstim.FepspSlope(key,sprintf('chan_num = %u',ch(i))),'fepsp_slope');
   tmp.slopes(i) = mean(slopes);
   tmp.slopes_med(i) = median(slopes);
   [tmp.freezing(i),tmp.shock(i)] = fetchn(beh.FearCond(key),'percent_freezing','shock_intensity');
end