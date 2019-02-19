function out =  check_for_files(excel_file, range, destination)
[~,~,d] = xlsread(excel_file, range);
d = d(:);
nFiles = length(d);
out = false(nFiles);

for i = 1:nFiles
    out(i) = exist(fullfile(destination, d{i}), 'dir');
    if ~out(i)
        disp(d{i})
    end
end



end