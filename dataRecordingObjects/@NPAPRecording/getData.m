function [V_uV,t_ms]=getData(obj,channels,startTime_ms,window_ms)

%obj,channels,startTime_ms,window_ms


%Extract recording data from file to memory
%Usage: [V_uV,t_ms]=obj.getData(channels,startTime_ms,window_ms);
%Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
%        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
%        window_ms - a scalar [1x1] with the window duration [ms].
%        name - the name of the recording (if empty takes the default name)
%Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
%        t_ms - A time vector relative to recording start (t=0 at start)


%%%%% IMPLEMENTATION:
binName = string(obj.dataFileNames(contains(obj.dataFileNames,'ap')));

path = obj.recordingDir;

meta = ReadMeta(binName, path);

nSamp = round((window_ms/1000)*str2double(meta.imSampRate));

V_uV = zeros(length(channels),length(startTime_ms),nSamp);

for trials = 1:length(startTime_ms)

    t0= startTime_ms(trials)/1000;
    samp0 = round(t0*str2double(meta.imSampRate));

    %[dataArray, channels] = ReadBin(samp0, nSamp, meta, binName,chanList, path);
    dataArray = ReadBin(samp0, nSamp, meta, binName, path);


    if obj.convertData2Double
        dataArray = GainCorrectIM(dataArray*1000000, channels, meta); %Multiply to convert to microvolts
        V_uV(:,trials,1:length(dataArray)) = dataArray(channels,:);
    else
        V_uV(:,trials,1:length(dataArray)) = dataArray(channels,:);
    end

    %convert to microvolts
end


t_ms = (1:nSamp)/obj.samplingFrequency*1000;

end

