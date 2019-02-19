function createSpikeDetectionSets(keys,detectMethodName)
% Populate detect.Params table for spike detection. This is the first step
% in setting up the spike detection process.
% MS 2016-07-28
for i = 1:length(keys)
    key = keys(i);
    dataDir = fetch1(acq.Ephys(key),'ephys_path');
    processedDir = strrep(dataDir,'raw','processed');
    key = fetch(acq.Ephys(key));
    key.detect_method_num = fetch1(detect.Methods(sprintf('detect_method_name = "%s"',...
        detectMethodName)),'detect_method_num');
    processedDir(1) = 'C';
    key.ephys_processed_path = processedDir;
    key.use_toolchain = 1;
    inserti(detect.Params,key)
end