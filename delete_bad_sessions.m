function delete_bad_sessions(sess_ts)

nSess = length(sess_ts);

for iSess = 1:nSess
    key = fetch(acq.Sessions(sprintf('session_path like "%%%s%%"',strtrim(sess_ts{iSess}))));
    x = acq.Sessions(key);
    x.del
end