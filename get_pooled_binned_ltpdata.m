function data = get_pooled_binned_ltpdata(ch_keys,datatype,args,varargin)
% function data = get_pooled_binned_ltpdata(ch_key,datatype,varargin)
% MS 2017-08-12
if isempty(ch_keys)
    data.t = [];
    data.val = [];
    return
end
tol = 5/60; % (in min) tolerange in bin center time jitter across sessions
r = struct;
nKeys = length(ch_keys);
for i = 1:nKeys
    ch_key = ch_keys(i);
    args.motion_idx_bound = args.all_motion_idx_bound(i,:);
    switch datatype
        case 'slopes'
            d = get_binned_epsp_resp(ch_key,args,varargin{:}); 
        case 'tdratio'
            d  = get_binned_tdrdata(ch_key,args,varargin{:});          
        case 'motion'
            d = get_binned_motion(ch_key,args,varargin{:});
        otherwise
            error('Given data type does not exist')
    end
    if ch_key.animal_id==175 && ch_key.session_start_time==5304402672
        d.bin_t(end) = d.bin_t(end-1)+diff(d.bin_t(1:2));
        warning('For mouse 175, the last bin was corrected')
    end
    r(i).val = d.bin_val';
    r(i).t = d.bin_t';
    r(i).len = length(d.bin_t);
end

% Pad with nan - pick timing from the longest vector
len = [r.len];
mdn = median(len);
% mdn = max(len);
if (mdn-round(mdn))~=0
    mdn = mode(len);
    disp(len)
    fprintf('median was %0.2f. So using mode value of %u instead\n',median(len),mdn)
end
ind = (len==mdn);

data.t = r(ind).t;
data.val = nan(mdn,nKeys);

for i = 1:nKeys
    sp = fetch1(acq.Sessions(ch_keys(i)),'session_path');
    aid = ch_keys(i).animal_id;
    rt = r(i).t;
    nT = length(rt);
    % Check each time bin for matching with the template time bin vector.
    % Then place the value in that matching location in the slope vector.
    % Allow some tolerance because across sessions, there may be a small
    % difference in the exact values of the bin center times.
    if mdn > nT
        for j = 1:nT
            ct = rt(j);
            [mv,idx] = min(abs(data.t-ct));
            if mv > tol
                fprintf('Bin center difference (%0.3f s) is above tolerance(%0.2f s) \n So the slope for time point %0.3f min was treated as NaN\n',mv*60,tol*60,ct)
                fprintf('Animal id: %u   session path: %s\n',aid,sp)
                val = nan;
            else
                val = r(i).val(j);
            end
            data.val(idx,i) = val;
        end
    else
        for j = 1:mdn
            ct = data.t(j);
            [mv,idx] = min(abs(rt-ct));
            if mv > tol
                fprintf('Bin center difference (%0.3f s) is above tolerance(%0.2f s) \n So the slope for time point %0.3f min was treated as NaN\n',mv*60,tol*60,ct)
                fprintf('Animal id: %u   session path: %s\n',aid,sp)
                val = nan;
            else
                val = r(i).val(idx);
            end
            data.val(j,i) = val;
        end
    end
end