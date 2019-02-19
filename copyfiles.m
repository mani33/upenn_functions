function copyfiles(Fp, source, destination)
% function copyfiles(Fp, source, destination)
% Inputs:
% Fp - session timestamp: eg.'2012-02-04_12-23-12'
% source - eg. W:\CheetahData
% destination - eg. y:\ephys\raw
%
% Vikas Chelur/Mani Subramaniyan

% Expected timestamp (Fp) pattern : eg. '2012-02-04_12-23-12'
expected = '[0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]{2}-[0-9]{2}-[0-9]{2}';
if ischar(Fp)
    Fp = {Fp};
end
source = strtrim(source);
destination = strtrim(destination);
nFiles = length(Fp);
if ~iscell(Fp)
    if ischar(Fp) && ~isempty(regexp(Fp,expected))%#ok
        Fp = strtrim(Fp);
        copyfile(fullfile(source, Fp), fullfile(destination, Fp));
    else
        disp('Empty folder supplied')
    end
else
    for i = 1:nFiles
        if ischar(Fp{i}) && ~isempty(regexp(Fp{i},expected))%#ok
            sess = strtrim(Fp{i});
            fsource = fullfile(source, sess);
            fdest = fullfile(destination, sess);
            if ~exist(fdest, 'dir')
                fprintf('Copying from %s to %s ... \n',fsource,fdest)
                copyfile(fsource, fdest);
                disp([sess '  successfully copied'])
            else
                list = dir(fsource);
                for j = 1:length(list)
                    ffdest = fullfile(destination,sess, list(j).name);
                    if ~exist(ffdest, 'file')
                        sf = fullfile(source, sess, list(j).name);
                        fprintf('Copying from %s to %s\n',sf,ffdest)
                        copyfile(sf, ffdest);
                    end
                end
                disp([sess '  already exists'])
            end
        else
            disp('Empty folder supplied')
        end
    end
end
