function [b_traces,bin_cen_t,trace_t] = getBinnedTraces( fpRespTraceRel, binwidth, varargin )
%getBinnedTraces get epsp traces averaged over a small time window
% 
% Inputs:
% key - database key that restricts the fepsp trace tuples
% chan - channel number
% binwidth - in minutes, to bin the responses
% Outputs:
% t - time in milli sec
% traces - raw averaged traces
% Mani Subramaniyan 2016-04-05
trace_bounds = [-5 25];

[et,y,t] = fetchn(fpRespTraceRel,'event_ts','y','t');


et = double(et);
[et,ind] = sort(et);
y = y(ind);
% All t's are pretty puch the same. so just pick one.
t = t{1}*1e-3; % convert to ms
tsInd = t >= trace_bounds(1) &  t <= trace_bounds(2);
trace_t = t(tsInd);
% Select the trace between given time points
tr = cellfun(@(y) y(tsInd), y, 'uni',false);
msz = min(cellfun(@length, tr));
tr = cellfun(@(x) x(1:msz), tr,'uni', false);
% Now, format time into minutes
et = format_time(et);
[b_traces,bin_cen_t] = bin_traces(tr,et,binwidth);


% subtract the baseline trace just after the stimulus artifact
% sel = (trace_t > 1.5) & (trace_t < 2.5);
sel = (trace_t > -4) & (trace_t < -1);
b_traces = cellfun(@(x) x-mean(x(sel)),b_traces,'uni',false);


