%example path:
% path = '\\132.66.45.127\data\Simon\Anesthesia_exp\Experiment_29_09_22\Insertion1\catgt_Experiment_29_09_22_1_g0';

function obj=extractMetaData(obj)

%Function to get NP metadata
path_r = obj.recordingDir;

fprintf('Extracting meta data from: %s...\n',path_r);

if iscell(path_r)
    path_r = cell2mat(path_r);
end


file = dir (path_r);
filenames = {file.name};
binFiles = filenames(contains(filenames,".bin"));

%binName, path

for stream = 1:length(binFiles)

    binName = string(binFiles(stream));

    meta = ReadMetaNP(binName,path_r);

    if contains(binName, '.ap.')
        metaAP = meta;
        metaAP.fileName = binName;
    elseif contains(binName, '.lf.')
        metaLF = meta;
        metaLF.fileName = binName;
    else
        metaNI = meta;
        metaNI.fileName = binName;
    end

end

clear meta

% Fill in properties
% 1. recordingName %(String) The name of the recording
out=regexp(metaAP.fileName,filesep,'split');
exp = erase(string(out(end)),'_g0_tcat.imec0.ap.bin');
obj.recordingName = char(exp);

%   %2. recordingDir % (String) Full directory containing the recorded session
%     out2=regexp(metaAP.fileName,"/"+exp,'split');
%     direc = string(out2(1));
%     obj.recordingDir = direc;

direc = path_r;
%3. dataFileNames % (Cell 1 x N)  array of N recording data file names
file = dir (direc);
filenames = {file.name};
binFiles = filenames(contains(filenames,".bin"));
obj.dataFileNames = binFiles;

%4. startDate %(1x1) Start date (time) of Recording (matlab date number format)
if isfield(metaAP, 'fileCreateTime_original')
    stTime = datetime(metaAP.fileCreateTime_original);
else
    stTime = datetime(metaAP.fileCreateTime);
end
obj.startDate = stTime;

%5. endDate %(1x1) End date (time) of Recording (matlab long format)

endTime = stTime + seconds(str2double(metaAP.fileTimeSecs));
obj.endDate = endTime;


%6. samplingFrequencyAP %(1xN) Sampling rate [Hz] of AP band
srateAP = str2double(metaAP.imSampRate);
obj.samplingFrequencyAP = srateAP;

%7. samplingFrequencyNI %(1xN) Sampling rate [Hz] of NIDQ
srateNI = str2double(metaNI.niSampRate);
obj.samplingFrequencyNI = srateNI;
obj.samplingFrequencyAnalog = srateNI;

%8. samplingFrequencyLF %(1xN) Sampling rate [Hz] of LF band
if exist("metaLF","var")
    srateLF = str2double(metaLF.imSampRate);
    obj.samplingFrequencyLF = srateLF;
else
    obj.samplingFrequencyLF = 0;
end

%6.0 samplingFrequency
obj.samplingFrequency = srateAP;

%9. recordingDuration_ms %(1x1) the total duration of the recording in [ms]
durMS = 1000*str2double(metaAP.fileTimeSecs);
obj.recordingDuration_ms = durMS;


%10. channelNumbers % (1xN) an array with integer channel numbers (>=1 integer)
if strcmp(metaAP.snsSaveChanSubset, 'all')
    chans = (1:str2double(metaAP.nSavedChans));
else
    chans = str2num(metaAP.snsSaveChanSubset);
    chans = chans + 1;
end
obj.channelNumbers = 1:numel(chans);

%10.0. channelNames % (Cell 1xN) a cell array with the N names of the channels

obj.channelNames= cellfun(@(x) char(string(x)),mat2cell(chans',ones(1,numel(chans))),'UniformOutput',0);

%triggerNames %the names of trigger channels (not critical)

%11. analogChannelNumbers % (1xN) the numbers of channels containing non-electrode analog inputs
if strcmp(metaNI.snsSaveChanSubset, 'all')
    chans1 = (1:str2double(metaNI.nSavedChans));
else
    chans1 = str2num(metaNI.snsSaveChanSubset);
    chans1 = chans1 + 1;
end
obj.analogChannelNumbers = chans1;

%12. analogChannelNames % (1xN) the numbers of channels containing non-electrode analog inputs

obj.analogChannelNames=cellfun(@(x) ['XA' char(string(x))],mat2cell(chans1(1:end)',ones(1,numel(chans1(1:end)))),'UniformOutput',0);
%chans1names(9) = "XD0";

%13.0

obj.chLayoutNumbers = metaAP.snsChanMap; %(MxN) The layout of the channel numbers in physical space arranged in an M by N grid

obj.chLayoutNames = (mat2cell(chans,1)); %(Cell MxN)The layout of the channel names in physical space arranged in an M by N grid

if isfield(metaAP,'imDatPrb_type')
    probeType = str2num(metaAP.imDatPrb_type);
else
    probeType = 0;
end


%14. layoutName %the name of the channel layout (electrode type)
if (probeType == 21) || (probeType == 24) || (probeType == 2013)
    obj.layoutName = 'NP_2.0_four_shanks';
else
    obj.layoutName = 'NP_1.0_staggered-3B1';
end

%13. ChannelMap

if isfield(metaAP,'snsGeomMap')
    [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = geomMapToGeom(metaAP);
elseif isfield(metaAP,'snsShankMap')
    [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = shankMapToGeom(metaAP);
end

nchans = numel(xCoord);
chanPos = zeros(2,nchans);
for i = 1:nchans
    currX = shankInd(i)*shankPitch + xCoord(i);
    %fprintf( fid, '%d\t%d\t%d\t%d\n', i-1, currX, yCoord(i), shankInd(i));

    chanPos(1,i) = currX;
    chanPos(2,i) = yCoord(i);
end

obj.chLayoutPositions = chanPos;

%Sync channel (bit) NI:

obj.syncNiChan = str2double(metaNI.syncNiChan);

%Digital channel NI:

obj.digChan = str2double(metaNI.nSavedChans)-1;

%calculate time conversion for syncing analog data and electrode data due to differences in sampling rates and drift between both traces [ms]
obj.syncSignalInElectrode = 1000*readmatrix([obj.recordingDir,filesep,dir(fullfile(obj.recordingDir, '*imec0.ap.xd_384_6_500.txt*')).name]);
obj.syncSignalInAnalog = 1000*readmatrix([obj.recordingDir,filesep,dir(fullfile(obj.recordingDir, sprintf('*nidq.xd_%d_%d_500.txt*',obj.digChan,obj.syncNiChan))).name]);


[APgain,LFgain] = ChanGainsIM(metaAP);
APgain= APgain(1); %takes gain of first channel. Asumes all channels have the same gain).

%15. MicrovoltsPerADimec % the digital to analog conversion value AP
obj.MicrovoltsPerAD = repmat(Int2Volts(metaAP)/APgain,1,length(chans)).*1000000;

%15.2 MicrovoltsPerADimec % the digital to analog conversion value LF
if isfield(metaAP,'imDatPrb_type')
    probeType = str2num(metaAP.imDatPrb_type);
else
    probeType = 0;
end

if (probeType == 21) || (probeType == 24) || (probeType == 2013)
    obj.MicrovoltsPerAD_LF=0;
else
    LFgain= LFgain(1);
    obj.MicrovoltsPerAD_LF = repmat(Int2Volts(metaAP)/LFgain,1,length(chans)).*1000000;
end

%16. MicrovoltsPerADnidq % the digital to analog conversion value
[MN,MA] = ChannelCountsNI(metaNI);
fI2V = Int2Volts(metaNI);

for i = 1:length(chans1)   % index into timepoint
    NIconv(i) = fI2V / ChanGainNI(i, MN, MA, metaNI);
end

convNidq = NIconv(1); %takes gain of first channel. Asumes all channels have the same gain).
obj.MicrovoltsPerADAnalog= convNidq.*1000000;



%17. Zero value for binary recordings
obj.ZeroADValue = 0;
obj.ZeroADValueAnalog = 0;

%18. Object type
obj.datatype = 'int16';        % class of data in the recording


%Channels
obj.nSavedChansNidq = metaNI.nSavedChans;
obj.nSavedChansImec = metaAP.nSavedChans;

obj.fileSizeBytesNI = metaNI.fileSizeBytes;
obj.fileSizeBytesAP = metaAP.fileSizeBytes;

if exist("metaLF","var")
    obj.fileSizeBytesLF = metaLF.fileSizeBytes;
else
    obj.fileSizeBytesLF = 0;
end

fprintf('saving meta data...\n');
obj.saveMetaData;


end % ReadMeta