function plot_all_chan_avg_stim_resp(key)
% Plots the average stimulation responses of all recorded channels.
% Example: plot_all_chan_avg_stim_resp('2017-07-24_11-07-29')

% Get key
% key = fetch(acq.Sessions(sprintf('session_path like "%%%s%%"',sess_timestamp)));
sess_timestamp = fetch1(acq.Sessions(key),'session_path');
if isempty(key)
    disp('No session found in database')
    return
end
% Get channels
ckeys = fetch(cont.Chan(key));
nChan = length(ckeys);
r = ceil(sqrt(nChan));

% Plot for each stimulating electrode separately (TTL 1 and 128)
stimTTL = unique(fetchn(acq.Events(ckeys,'event_ttl>0'),'event_ttl'));
nStim = length(stimTTL);

for iStim = 1:nStim
    ttl = stimTTL(iStim);
    figure
    h = nan(1,1);
    hi = 0;
    for iChan = 1:nChan
        ckey = ckeys(iChan);
        hi = hi+1;
        [chname,cn] = fetchn(cont.Chan(ckey),'chan_name','chan_num');
        
        [y,t] = fetchn(cstim.FpRespTrace(ckey,sprintf('event_ttl = %u',ttl)),'y','t');
        if isempty(y)
            fprintf('For chan %u, FpRespTrace tuples don''nt exist\n',cn)
        else
            h(hi) = subplot(r,r,iChan); box off
        minLen = min(cellfun(@length,y));
        ytrim = cellfun(@(x) x(1:minLen),y,'uni',false);
        ytrim = [ytrim{:}];
        t = t{1};
        t = t(1:minLen)/1000; % msec
        yavg = mean(ytrim,2)*1000; % mV
        plot(t,yavg);
        axis tight; box off
        
        title(sprintf('%s: %u',chname{:},cn))
        end
    end
    linkaxes(h,'x')
    suptitle(['mouse: ' num2str(ckeys(iChan).animal_id) ' ' sess_timestamp sprintf('  : TTL %u', ttl)])
end
