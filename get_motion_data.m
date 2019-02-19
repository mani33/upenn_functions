function [t, m] = get_motion_data(key)
assert(length(key)==1,'Only one key at a time please')
mobj = cont.Motion(key);
if count(mobj)==0
    mobj = cont.MotionInVid(key);
    if count(mobj)==0
        disp('No motion data found')
        t = nan;
        m = nan;
        return
    end
end
% cont.Motion object already trims the motion data to bracket the first
% and last electrical pulse.
[t,m] = fetchn(mobj,'t','dist_var'); % time in sec

t = t{:}'; % in us
t = (t-t(1))*(1e-6)/60; % in min
m = m{:}';
