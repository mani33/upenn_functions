function copy_data(excelFile,range,data_source_path,data_dest_path)
% Copies data from data source to data destination. 
% excelFile - file that contains just the timestamp of the recording
% folder.

[~,~,d] = xlsread(excelFile,1, range);
d = d(:);
copyfiles(d, data_source_path,data_dest_path)
