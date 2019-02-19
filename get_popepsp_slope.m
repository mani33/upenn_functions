function s = get_popepsp_slope(t,y,varargin)

% Use only microseconds as time unit except for convolution
args.bounds = [3 15]*1000;% 
args.slope_win = 1*1000; % slope window
args.std = 0.5; % gauss window std in ms
args = parseVarArgs(args,varargin{:});

dt = diff(t(1:2));
gw = getGausswin(args.std,dt/1000);

% smooth the data to get the slopes without noise
yso = mconv(y,gw);
% Compute slope
dy = diff(yso);
% Clip the slope
sel = (t>=args.bounds(1)) & t<args.bounds(2);
dys = dy;
dys(~sel) = nan;

% In this piece, find the window to calculate slope
[~,ind] = max(abs(dys+1e12));
ns = round(args.slope_win/dt);
sind = ind + (round(-ns/2):round(ns/2));
ys = y(sind);
ts = t(sind); % in seconds

% Fit slope
X = [ts,ones(size(ts))];
B = regress(ys,X);
s = B(1);
yi = X*B;
% plot(X(:,1)*1e6,yi,'m*')
% plot([X(1,1) X(1,1)]*1e6,ylim,'k--')
% plot([X(end,1) X(end,1)]*1e6,ylim,'k--')






          