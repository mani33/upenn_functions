function [m_data,tm,td_data,t_td] = get_mouse_motion_and_tdratio_for_ltp(chkeys_pre,chkeys_post,varargin)
% function get_mouse_motion_and_tdratio_for_ltp(keys_pre,keys_post)
% Get the mouse movement displacement and theta_delta ratio for baseline
% and post-treatment ltp experiment. If more than one pair of keys
% provided, averaging will be done, optionally. Keys should have ch num,
% for theta delta ratio is available for individual channels.
%
% example: for single mouse
% get_mouse_motion_and_tdratio_for_ltp(chkey_pre,chkey_post)
% for averaging across 4 mice
% get_mouse_motion_and_tdratio_for_ltp(chkeys_pre,chkeys_post,'avg',1)
%
% MS 2017-08-03
%
args.avg = 1;

args = parseVarArgs(args,varargin{:});

tmp = struct;
n = length(chkeys_pre);

for i = 1:n
    key_pre = chkeys_pre(i);
    key_post = chkeys_post(i);
    
    %     [tb_pre,te_pre] = get_trimming_times(key_pre);
    %     [tb_post,te_post] = get_trimming_times(key_post);
    
    [t, tmp.m_pre{i}] = get_motion_data(key_pre);
    tmp.tm_pre{i} = t-t(end);
    [tmp.tm_post{i}, tmp.m_post{i}] = get_motion_data(key_post);
    
    [t, tmp.td_pre{i}] = get_tdratio_data(key_pre);
    tmp.tt_pre{i} = t - t(end);
    [tmp.tt_post{i}, tmp.td_post{i}] = get_tdratio_data(key_post);
    
end

% average if necessary

tlen_pre = cellfun(@length, tmp.tm_pre);
tlen_post = cellfun(@length, tmp.tm_post);
md1 = median(tlen_pre);
md2 = median(tlen_post);

% pad times/resp for pre
padsz = arrayfun(@(x) x, md1-tlen_pre,'uni',0);
tm_pre = cellfun(@(x,np) [nan(np,1); x],tmp.tm_pre,padsz,'uni',0);
mpre = cellfun(@(x,np) [nan(np,1); x],tmp.m_pre,padsz,'uni',0);

% Trim post resp since one session had more than 180 pulses
tm_post = cellfun(@(x) x(1:md2),tmp.tm_post,'uni',0);
mpost = cellfun(@(x) x(1:md2),tmp.m_post,'uni',0);

% td ratio
tdlen_pre = cellfun(@length, tmp.tt_pre);
tdlen_post = cellfun(@length, tmp.tt_post);
md1 = median(tdlen_pre);
md2 = median(tdlen_post);

% pad times/resp for pre
padsz = arrayfun(@(x) x, md1-tdlen_pre,'uni',0);
tt_pre = cellfun(@(x,np) [nan(np,1); x],tmp.tt_pre,padsz,'uni',0);
tdpre = cellfun(@(x,np) [nan(np,1); x],tmp.td_pre,padsz,'uni',0);

% Trim post resp since one session had more than 180 pulses
tt_post = cellfun(@(x) x(1:md2),tmp.tt_post,'uni',0);
tdpost = cellfun(@(x) x(1:md2),tmp.td_post,'uni',0);

% Average now
m_data = [[mpre{:}];[mpost{:}]];
tm = [[tm_pre{:}];[tm_post{:}]];
td_data = [[tdpre{:}];[tdpost{:}]];
t_td = [[tt_pre{:}];[tt_post{:}]];


%     m_data = [nanmean([mpre{:}],2); nanmean([mpost{:}],2)];
%     tm = [nanmean([tm_pre{:}],2); nanmean([tm_post{:}],2)];
%     td_data = [nanmean([tdpre{:}],2); nanmean([tdpost{:}],2)];
%     t_td = [nanmean([tt_pre{:}],2); nanmean([tt_post{:}],2)];

if args.avg && n > 1
    m_data = nanmean(m_data,2);
    tm = nanmean(tm,2);
    td_data = nanmean(td_data,2);
    t_td = nanmean(t_td,2);
end





% function [tb,te] = get_trimming_times(key)
% % Since motion tracking starts before the first electrical pulse and stops
% % after the last pulse, we will need to trim the motion data.
% ets = double(fetchn(acq.Events(key,'event_ttl in (1,128)'),'event_ts'));
% epd = fetch(acq.Ephys(key),'*');
% tb = (min(ets) - double(epd.ephys_start_time))*(1e-6)/60; % min
% stim_dur = (max(ets)-min(ets))*(1e-6)/60; % min
% te = tb + stim_dur;