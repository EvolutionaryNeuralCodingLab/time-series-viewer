function [V_uV,t_ms]=getAnalogData(obj,channels,startTime_ms,window_ms)
%Extract Neuralynx recording data from file to memory
%Usage: [V_uV,t_ms]=obj.getData(channels,startTime_ms,window_ms);
%Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
%        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
%        window_ms - a scalar [1x1] with the window duration [ms].
%Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
%        t_ms - A time vector relative to recording start (t=0 at start)

%transpose startTime_ms if it was provided a column vector and not as row vector (to prevent an error further on)
if size(startTime_ms,1)>1
    startTime_ms=startTime_ms';
end

%Implementation:
% Todo: Place all the opertations below in the class constructor and run only once
binName = string(obj.dataFileNames(contains(obj.dataFileNames,'nidq')));

path = obj.recordingDir;

meta = ReadMeta(binName, path);

chanList = channels;

nChan = str2double(meta.nSavedChans);

nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);

if obj.syncAnalogWithElectrode
%determine start and end times in analog recording based on electrode recordings using the sync data in syncSignal
    [start_End] = interp1(obj.syncSignalInElectrode, obj.syncSignalInAnalog, [startTime_ms; startTime_ms+window_ms], 'linear');
else
    start_End=[startTime_ms;startTime_ms+window_ms];
end
nSamp = round(((start_End(2,1)-start_End(1,1))/1000)*str2double(meta.niSampRate));
nSamp = min(nSamp, nFileSamp);


V_uV = zeros(numel(channels),numel(startTime_ms),nSamp);

for trials = 1:length(startTime_ms)

    t0= start_End(1,trials)/1000;

    samp0 = round(t0*str2double(meta.niSampRate));
% Todo: ReadBin reads all the channels - should be modfied to read only the relevant channels (at least in the simple case of 1 channel)
    dataArray = ReadBin(samp0, nSamp, meta, binName, path);

    dataArray = dataArray(chanList,:);

    if obj.convertData2Double
        V_uV(:,trials,:) = GainCorrectNI(dataArray, chanList, meta); %convert to microvolts
    else
        V_uV = dataArray;
    end

end

%modify time vector based on the sync signal lag of the first trial
t_ms = (1:size(V_uV,3))/obj.samplingFrequencyNI*1000*((start_End(2,1)-start_End(1,1))/window_ms);


%Added by Mark 28/6/25 - the code below resulted in size differences between V_uV and t_ms
%t_ms = 0:str2double(meta.niSampRate)/1000:str2double(meta.fileTimeSecs)*1000;
%t_ms= round(t_ms*1e6)/1e6;
%digits(6)

end
