function [p,yhat,rsq] = test_for_baseline_stability(x,y)
% function [p,yhat,rsq] = test_for_baseline_stability(x,y)
% Tests if the baseline data points are flat. 
% Inputs:
% x - independent variable (time)
% y - dependent variable (slope)
% Mani Subramaniyan: 2018-01-18

lm = fitlm(x(:),y(:),'linear','RobustOpts','off');
yhat = predict(lm,x(:));
p = coefTest(lm);
rsq = roundn(lm.Rsquared.Ordinary,-2);

