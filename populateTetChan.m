function populateTetChan(keys,varargin)
% populateTetChan(keys,param1,paramVal1,param2,paramVal2,...)
%-----------------------------------------------------------------------------------------
% POPULATETETCHAN - populate cont.TetChan table
%
% example: populateTetChan(keys)
%
% This function is called by:
% This function calls:
% MAT-files required:
%
% See also:

% Author: Mani Subramaniyan
% Date created: 2015-10-25
% Last revision: 2015-10-25
% Created in Matlab version: 8.3.0.532 (R2014a)
%-----------------------------------------------------------------------------------------

% We assume the following mapping
nTet = 4;

% Find the channel number and assign tetrode number based on the above
% mapping
keys = fetch(acq.Ephys(keys));
nk = length(keys);
for i = 1:nk
    key = keys(i);
    ck = fetch(cont.Chan(key));
    cn = [ck.chan_num];
    for tt = 1:nTet
        ch = (0:3)+(tt-1)*4;
        [ic,ind] = ismember(ch,cn);
        assert(length(find(ic))==4,'All 4 channels not available in cont.Chan table for given tetrode')
        key.tet_num = tt;
        key.chan_nums = sort(cn(ind));
        inserti(cont.TetChan,key)
    end
end








