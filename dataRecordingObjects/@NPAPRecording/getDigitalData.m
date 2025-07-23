function [Dig,T_ms]=getDigitalData(obj,startTime_ms,window_ms)
%Extract MCRack digital data from file to memory
%Usage: [Dig,T_ms]=getDigitalData(startTime_ms,window_ms,name)
%Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
%        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
%        window_ms - a scalar [1x1] with the window duration [ms].
%        name - the name of the stream (if not entered, default name is used)
%Output: Dig - A 3D matrix [nChannels x nTrials x nSamples] with digitalData waveforms across specified channels and trials
%        T_ms - A time vector relative to recording start (t=0 at start)


binName = string(obj.dataFileNames(contains(obj.dataFileNames,'nidq')));

path = obj.recordingDir;

meta = ReadMeta(binName, path);

dw = 8;
dLineList = [1, 2, 3];

channels =1;

nSamp = round((window_ms/1000)*str2double(meta.imSampRate));

Dig = zeros(length(channels),length(startTime_ms),nSamp);

for trials = 1:length(startTime_ms)

    t0= startTime_ms(trials)/1000;
    samp0 = round(t0*str2double(meta.imSampRate));

    %[dataArray, channels] = ReadBin(samp0, nSamp, meta, binName,chanList, path);
    dataArray = ReadBin(samp0, nSamp, meta, binName, path);



    dataArray = GainCorrectNI(dataArray*1000000, channels, meta); %Multiply to convert to microvolts
    Dig(:,trials,1:length(dataArray)) = dataArray(channels,:);

    %convert to microvolts
end


T_ms = (1:nSamp)/obj.samplingFrequency*1000;

end

