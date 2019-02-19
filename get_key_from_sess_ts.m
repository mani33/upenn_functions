function keys = get_key_from_sess_ts(sess_timestamps)

if ischar(sess_timestamps)
    sess_timestamps = {sess_timestamps};
end

sess_timestamps = sess_timestamps(:);
n = length(sess_timestamps);
k = 0;

for i = 1:n
    st = sess_timestamps{i};
    if ~isnan(st)
        ckey = fetch(acq.Sessions(sprintf('session_path like "%%%s%%"',st)));
        if ~isempty(ckey.animal_id)
            k = k + 1;
            keys(k) = ckey;
        else
            fprintf('No key found for session: %s\n',st)
        end
    end
end