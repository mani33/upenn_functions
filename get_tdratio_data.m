function [t, td] = get_tdratio_data(key)

mobj = cont.TDratio(key);
if count(mobj)~=0
    % cont.TDratio object already trims the motion data to bracket the first
    % and last electrical pulse.
    [t,td] = fetchn(mobj,'t','td_ratio'); % time in us
    t = t{:}';
    t = t-t(1);
    t = t*(1e-6)/60; % in min
    td = td{:}';
else
    error('No td ratio data found')
end