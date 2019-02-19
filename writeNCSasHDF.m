function writeNCSasHDF(sourceFile,outFile)
% Rewrite Neuralynx continuous data files(*.ncs)file in hdf5 format.
% MS 2015-04-22
% Load data

br = baseReaderNeuralynx(sourceFile);
Fs = getSamplingRate(br);
rs = getRecordSize(br);
n = getNbSamples(br);
nRecords = round(n/rs);
chName = getChannelNames(br);
scale = getScale(br);



% Create a local file because writing over the network causes problems
if strfind(outFile,'y:') || strfind(outFile,'Y:')
    outFile = strrep(outFile,'y:','C:');
    mkdir(fileparts(outFile))
else
    error('The out file location is not y or Y drive')
end

fp = H5Tools.createFile(outFile, 'driver', 'family');


for iRec = 1:nRecords
    [timeStamps,nValidSamples,records] =  Nlx2MatCSC(sourceFile,[1 0 0 1 1],0,2,iRec); % recordStartTimes
   % scale to volts
    order = numel(scale);
    if order == 1
        x = x * br.scale;
    else
        y = 0;
        for i = 1:order
            y = y + x.^(i - 1) * br.scale(i);
        end
        x = y;
    end
    
    % write data to disc
    if iRec == 1
        [dataSet, written] = seedDataset(fp,x);
    else
        written = written + extendDataset(dataSet, x, written);
    end
    
    displayProgress(iRec,nRecords)
end

H5D.close(dataSet);


% Now create/copy the remaining attributes etc.
H5Tools.writeAttribute(fp, 'Fs', Fs);
H5Tools.writeAttribute(fp, 'channelNames', chName);
H5Tools.writeAttribute(fp, 'class', 'Electrophysiology');
H5Tools.writeAttribute(fp, 'version', 1);
H5Tools.writeAttribute(fp, 'scale', 1);  % data are in Volts
H5Tools.writeAttribute(fp, 't0', getT0(br));
H5F.close(fp);



function [dataSet, written] = seedDataset(fp, data)

nbDims = 2;
dataDims = size(data);
dataDims(1:2) = dataDims([2 1]);
dataType = H5Tools.getHDF5Type(data);
dataSpace = H5S.create_simple(nbDims, dataDims, [dataDims(1) -1]);

setProps = H5P.create('H5P_DATASET_CREATE'); % create property list
chunkSize = [4, 100000]; 		% define chunk size
chunkSize = min(chunkSize, dataDims);
H5P.set_chunk(setProps, chunkSize); % set chunk size in property list

dataSet = H5D.create(fp, '/data', dataType, dataSpace, setProps);
H5D.write(dataSet, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', data);

H5P.close(setProps);
H5T.close(dataType);
H5S.close(dataSpace);
written = size(data, 1);


function written = extendDataset(dataSet, data, written)

% Extend dataset
H5D.extend(dataSet, [size(data,2), written+size(data,1)])

% Select appended part of the dataset
fileSpace = H5D.get_space(dataSet);
H5S.select_hyperslab(fileSpace, 'H5S_SELECT_SET', [0, written], [], fliplr(size(data)), []);

% Create a memory dataspace of equal size.
memSpace = H5S.create_simple(2, fliplr(size(data)), []);

% And write the data
H5D.write(dataSet, 'H5ML_DEFAULT', memSpace, fileSpace, 'H5P_DEFAULT', data);

% Clean up
H5S.close(memSpace);
H5S.close(fileSpace);
written = size(data,1);



function decFactors = calcDecimationFactors(samplingRate, cutoffFreq)

targetRate = cutoffFreq * 2 / 0.8;		% 0.8 is the cutoff of the Chebyshev filter in decimate()
coeff = samplingRate / targetRate;

if (coeff < 2)
    error('Cannot decimate by %g', coeff)
end

% Calculate a series of decimation factors
decFactors = [];
testFactors = 13:-1:2;
while (coeff > 13)
    rems = mod(coeff, testFactors);
    [~, ix] = min(rems);
    decFactors = [decFactors, testFactors(ix)]; %#ok<AGROW>
    coeff = coeff / testFactors(ix);
end

coeff = floor(coeff);
if (coeff >= 2)
    decFactors = [decFactors, coeff];
end


function out = decimatePackage(data, factor)

[m,n] = size(data);
% crop package at a multiple of decimation factor. this is important
% because otherwise decimate will cause random jitter of up to one sample
m = fix(m / factor) * factor; 
out = zeros(m/factor,n);
for col = 1:n
    out(:,col) = decimate(data(1:m,col), factor);
end
