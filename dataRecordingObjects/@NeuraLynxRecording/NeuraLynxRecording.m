classdef NeuraLynxRecording < dataRecording
    properties
        dataFileNumbers %Numbers of recorded channel files
        triggerFileName %The name (including full path) of the trigger file in Neuralynx recording
        firstRecordStartTime_us %The time stamp of the first recorded waveform
        lastRecordStartTime_us %The time stamp of the last recorded waveform
        nRecordsPerChannel %the number of records recorded for every channel
        nonActiveChannelNumbers %channel numbers with no records (turned off during experiment)
        A2DMaxValue %
        A2DBitVolts %Quantization of data in terms of voltage (coversion between bits and volts)
        saturationVoltage_uV %Max voltage allowed by the system
        recordLength_us
        sample_us
        nRecords
        maxTTLBit = 9;
        invertTraces = 1;
        %{
        Properties of parent class: dataRecording
        recordingDir %Full directory containing the recorded session
        startDate %Start date (time) of Recording (matlab long format)
        endDate %End date (time) of Recording (matlab long format)
        channelNames %a cell array with the names of the channels
        channelNumbers %an array with integer channel numbers
        recordingDuration_ms %The end time of the recording (including all 512 samples of the last record)
        samplingFrequency %Sampling rate
        dataFileNames %Name of recording data files
        dspLowCutFrequency %Low-pass cutoff frequency in the Neuralynx DSP (in raw data)
        dspHighCutFrequency %High-pass cutoff frequency in the Neuralynx DSP (in raw data)
        %}
    end
    properties (Constant)
        DefaultBitsForVisualStimTriggers=[1]; %Default TTL bits extracted
        defaultChannelNumberForHearderExtraction=1; %Channel from which meta data is extracted (assumes that all channel were recorded similarly
        samplesPerRecord=512; %Number of samples in one record of neuralynx data
        listOfFieldsToExtractFromHeader=...
            {'RecordSize','SamplingFrequency','ADMaxValue','ADBitVolts','DspLowCutFrequency','DspHighCutFrequency';
            'nRecords','samplingFrequency','A2DMaxValue','A2DBitVolts','dspLowCutFrequency','dspHighCutFrequency'}; %List of fields to be extracted from the header and their associated property names
        defaultLocalDir='C:\Users\Tulip\Documents\Academic\Post-Doc\Experiments'; %Default directory from which search starts
        fileExtension='ncs';
        numberOfCharFromEndToBaseName = 7;
        signalBits = 16;
    end
    methods
        function obj=reloadMetaData(obj)
            %Reloads the meta data from directory
            %Usage: obj.reloadMetaData();
            %Input: []
            obj=getChannelInformation(obj);
            obj=getHeaderInformation(obj);
        end
        
        
        
        function [V_uV,t_ms]=getData(obj,channels,startTime_ms,window_ms)
            %Extract Neuralynx recording data from file to memory
            %Usage: [V_uV,t_ms]=obj.getData(channels,startTime_ms,window_ms);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        t_ms - A time vector relative to recording start (t=0 at start)
            startTime_us=startTime_ms*1e3+obj.firstRecordStartTime_us;
            windowSamples=round(window_ms/1e3*obj.samplingFrequency);
            window_us=window_ms*1e3;
            nWindows=numel(startTime_us);
            
            if isempty(channels) %if no channels are entered, get all channels
                channels=obj.channelNumbers;
            end
            nCh=numel(channels);
            
            V_uV=NaN(windowSamples,nWindows,nCh); %initialize waveform matrix
            
            %generate time stamps for block waveform extraction and extract waveforms from file
            singleTrialTimeStamps=(0:obj.recordLength_us:(window_us+obj.recordLength_us-obj.sample_us))';
            timeVector=bsxfun(@plus,singleTrialTimeStamps,startTime_us);
            
                  
            pStampOutsideRange=timeVector<obj.firstRecordStartTime_us | timeVector>obj.lastRecordStartTime_us-obj.recordLength_us;
            completeTrials=~any(pStampOutsideRange);
            %go over trials that are not cut
            
            if any(completeTrials) %check for existance of full trials
                startTimeFull_us=startTime_us(completeTrials);
                timeVectorFull=timeVector(:,completeTrials);
                startOffset=(find(completeTrials,1,'first')-1)*windowSamples;
                nWindowsFull=sum(completeTrials);
                
                nSingleTrialTimeStamps=numel(singleTrialTimeStamps);
                startTimeRecords_uS=bsxfun(@times,ones(nSingleTrialTimeStamps,1),startTimeFull_us);
                
                %run over all channels
                for i=1:nCh
                    %extract Neuralynx data
                    [timeStampsData, dataSamples] = Nlx2MatCSC([obj.recordingDir '\' obj.dataFileNames{channels(i)}],[1 0 0 0 1],0,5,timeVectorFull(:));
                    %add functionality to NlxCSC2mat for extracting mode 5
                    
                    try
                        idx=round(bsxfun(@plus,(timeStampsData-startTimeRecords_uS(:)')/obj.sample_us,(1:obj.samplesPerRecord)'));%indices in V_uV but including samples outside range
                        pIdx=find(idx>0 & idx<=windowSamples); %find the indices that will be within the limits in V_uV
                        cumSample=((0:windowSamples:(windowSamples*(nWindowsFull)-1))'*ones(1,nSingleTrialTimeStamps))';
                        idx=bsxfun(@plus,idx,cumSample(:)');
                        
                        %convert to [Volt] units and generate time vector [ms]
                        V_uV(idx(pIdx)+startOffset+(i-1)*windowSamples*nWindows)=dataSamples(pIdx);
                        
                    catch  %weird bug where it cant find timestamps in mode 5, but there is not any missing...
                        startTimeIndex_us=startTime_ms*1e3;
                        
                        singleTrialTimeStampNumber=(0:floor(window_us/obj.recordLength_us)+1)';
                        indexVector=bsxfun(@plus,singleTrialTimeStampNumber,floor((startTime_us-startTime_us(1))/obj.recordLength_us)+1);
                        
                        indexVectorFull=indexVector(:,completeTrials);
                        [timeStampsD, dataSamples] = Nlx2MatCSC([obj.recordingDir '\' obj.dataFileNames{channels(i)}],[1 0 0 0 1],0,4,[startTimeFull_us(1) (startTimeFull_us(end)+window_us)]);
                        timeStampsData=timeStampsD(indexVectorFull(:));
                        dataSamples=dataSamples(:,indexVectorFull(:));
                        idx=round(bsxfun(@plus,(timeStampsData-startTimeRecords_uS(:)')/obj.sample_us,(1:obj.samplesPerRecord)'));%indices in V_uV but including samples outside range
                        pIdx=find(idx>0 & idx<=windowSamples); %find the indices that will be within the limits in V_uV
                        cumSample=((0:windowSamples:(windowSamples*(nWindowsFull)-1))'*ones(1,nSingleTrialTimeStamps))';
                        idx=bsxfun(@plus,idx,cumSample(:)');
                        %convert to [Volt] units and generate time vector [ms]
                        V_uV(idx(pIdx)+startOffset+(i-1)*windowSamples*nWindowsFull)=dataSamples(pIdx);
                    end
                end
            end
            
            if any(~completeTrials)
                
                for j=find(~completeTrials)
                    
                    timeVectorCut=timeVector(~pStampOutsideRange(:,j),j);
                    nSingleTrialTimeStamps=numel(timeVectorCut);   
                  windowSamplesCut=length(timeVectorCut)*obj.samplesPerRecord;
                   V_uV=NaN(length(timeVectorCut)*obj.samplesPerRecord,nWindows,nCh); %initialize waveform matrix with cut times
                 
                    
                    startTimeRecords_uS=bsxfun(@times,ones(nSingleTrialTimeStamps,1),startTime_us(j));
                    %check if some time stamps exceed recording duration
                   
                    %run over all channels
                    for i=1:nCh
                        %extract Neuralynx data
                        [timeStampsData, dataSamples] = Nlx2MatCSC([obj.recordingDir '\' obj.dataFileNames{channels(i)}],[1 0 0 0 1],0,5,timeVectorCut(:));
                        
                        idx=round(bsxfun(@plus,(timeStampsData-startTimeRecords_uS(:)')/obj.sample_us,(1:obj.samplesPerRecord)'));%indices in V_uV but including samples outside range
                        pIdx=find(idx>0 & idx<=windowSamples); %find the indices that will be within the limits in V_uV
                        
                        cumSample=((0:windowSamplesCut:(windowSamplesCut-1))'*ones(1,nSingleTrialTimeStamps))';
                        idx=bsxfun(@plus,idx,cumSample(:)');
                        
                        %convert to [Volt] units and generate time vector [ms]
                        V_uV(idx(pIdx)  +  (j-1)* windowSamplesCut +  (i-1)* windowSamplesCut*nWindows)=dataSamples(pIdx);
                    end
                end
            end
            
            V_uV=(2*(0.5-obj.invertTraces))*permute(V_uV,[3 2 1])*(obj.A2DBitVolts*1e6); %first part controls the signal inversion
      
         
            if nargout==2
                t_ms=(1:windowSamples)*(1e3/obj.samplingFrequency);
            end
        end
        
        function [T_ms,chNum,chName]=getTrigger(obj,startTime_ms,window_ms,name)
            %Extract triggers from file Neuralynx recording
            %Usage : [T_ms]=obj.getTrigger(name,startTime_ms,endTime_ms)
            %Input : name - which bit to extract for time stamps (out of 8,default = first bit, 1)
            %        startTime_ms - start time [ms].
            %        window_ms - the window duration [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %
            %Output: T_ms - trigger times [ms] - different triggers are arranged in a cell array
            
            %read Data from *.nev files (NeuroLynx) - files containing the data
            %fieldSelectionTrigger: (1): Timestamps, (2): Event IDs, (3): TTLs, (4): Extras, (5): Event Strings
            %Bits are coded as decimal integers
            if nargin>=3
                [timeStampsTrig, TTLsTrig,eventString] = Nlx2MatEV...
                    ([obj.recordingDir '\' obj.triggerFileName],[1 0 1 0 1],0,4,1e3*[startTime_ms; startTime_ms+window_ms]+obj.firstRecordStartTime_us);
            else
                [timeStampsTrig, TTLsTrig,eventString] = Nlx2MatEV...
                    ([obj.recordingDir '\' obj.triggerFileName],[1 0 1 0 1],0,1,[]);
            end
            
            %calculate bit from integer sequence
            maxTTLBit=8;
            %[~,e]=log2(max(TTLsTrig));
            pSimpleTTLBinary=TTLsTrig<=2^(maxTTLBit);
            pSimpleTTL=find(pSimpleTTLBinary);
            TTLBinaryCh=logical(fliplr(rem(floor(TTLsTrig(pSimpleTTL)'*pow2(1-maxTTLBit:0)),2)));
            
            %check if specific bits were selected for extraction
            if nargin~=4
                bits2Check=1:maxTTLBit;
            else
                bits2Check=name;
            end
            
            %extract bits time stamps and names
            for i=bits2Check%(sum(TTLBinaryCh)>=1)
                pTmp=pSimpleTTL(1+find(diff(TTLBinaryCh(:,i))~=0));
                T_ms{i}=(   timeStampsTrig(pTmp) - obj.firstRecordStartTime_us   )/1e3;
                chName{i}=eventString(pTmp);
                chNumber(i)=i;
            end
            
            %extract TTLs that were coded as integers
            if isempty(bits2Check)
                bits2Check=0;
            end
            complexTTLStart=max(bits2Check)+1;
            complexTTLs=unique(TTLsTrig(~pSimpleTTLBinary));
            for i=1:numel(complexTTLs)
                pTmp=(find(TTLsTrig==complexTTLs(i)));
                T_ms{i+complexTTLStart}=(timeStampsTrig(pTmp)-obj.firstRecordStartTime_us)/1e3;
                chName{i+complexTTLStart}=eventString(pTmp);
                chNumber(i+complexTTLStart)=complexTTLs(i);
            end
            
        end
    end
    
    methods (Hidden)
        %class constructor
        function obj = NeuraLynxRecording(recordingFile)
            %Usage: obj = NeuraLynxRecording(recordingFile),
            
            %get data files
            if nargin==0
                recordingFile=[];
            elseif nargin>1
                disp('Object was not constructed since too many parameters were given at construction');
                return;
            end
            obj.folderMode=true;
            obj=obj.getRecordingFiles(recordingFile);
            obj.recordingDir=[obj.recordingDir filesep obj.recordingName];
            
            metaDataFile=[obj.recordingDir filesep 'metaData.mat'];
            if exist(metaDataFile,'file')
                obj=loadMetaData(obj,metaDataFile);
            else
                obj=getChannelInformation(obj);
                obj=getHeaderInformation(obj);
                triggerNames=dir([obj.recordingDir filesep '*.nev']);
                if numel(triggerNames)~=1
                    disp('triggers not extracted - more than one trigger file in directory');
                else
                    obj.triggerFileName=triggerNames.name;
                end
                obj.saveMetaData;
            end
        end
        
        function obj=getChannelInformation(obj)
            %Reads recording files from recording directory and extracts file information
            allFiles=dir([obj.recordingDir filesep '*.ncs']);
            unsorteddataFileNames={allFiles.name};
            unsorteddataFileNumbers=arrayfun(@(x) str2double(x.name(4:end-4)), allFiles);
            
            [~,places]=sort(unsorteddataFileNumbers);
            obj.dataFileNames=unsorteddataFileNames(places);
            obj.dataFileNumbers=unsorteddataFileNumbers(places);
            obj.channelNumbers=obj.dataFileNumbers;
            h = waitbar(0,'Initializing NeuraLynx object data (performed only on first access)...');
            nCh=numel(obj.channelNumbers);
            for i=1:nCh
                waitbar(i / nCh)
                try
                    records=Nlx2MatCSC_v3([obj.recordingDir filesep obj.dataFileNames{i}],[0 0 0 1 0],0,1,[]);
                catch
                    records=[];
                end
                %records=Mat2NlxCSC([obj.recordingDir filesep obj.dataFileNames{i}],[0 0 0 1 0],0,1,[]);
                %Check if using Nlx2MatCSC_v3 is faster
                obj.nRecordsPerChannel(i) = numel(records);
            end
            close(h);
            p=obj.nRecordsPerChannel>1;
            obj.nonActiveChannelNumbers=obj.channelNumbers(~p);
            obj.channelNumbers=obj.channelNumbers(p);
            obj.channelNames=num2cell(obj.channelNumbers);
            obj.nRecordsPerChannel=obj.nRecordsPerChannel(p);
            if ~all(obj.nRecordsPerChannel==obj.nRecordsPerChannel(1))
                warning('Notice that a different number of records was collected for different channels');
            end
        end
        
        function obj=getHeaderInformation(obj) %get parameters from file header
            %get header and first time stamp from recordings
            [timeStamps, header] = Nlx2MatCSC...
                ([obj.recordingDir filesep obj.dataFileNames{obj.defaultChannelNumberForHearderExtraction}],[1 0 0 0 0],1,1,[]);
            
            for i=1:numel(header)
                if ~isempty(header{i})
                    if header{i}(1)=='-'
                        pValue=find(header{i}==' ',1,'first');
                        paramLineInHeader=find(strcmp(header{i}(2:pValue-1),obj.listOfFieldsToExtractFromHeader(1,:)));
                        if ~isempty(paramLineInHeader)
                            obj.(obj.listOfFieldsToExtractFromHeader{2,paramLineInHeader})=...
                                str2num(header{i}(pValue+1:end));
                        end
                    elseif header{i}(1:14)=='## Time Opened'
                        pValue=find(header{i}==' ',7,'first');
                        obj.startDate=datenum(header{i}([pValue(4)+1:pValue(5) pValue(7)+1:end]));
                    elseif header{i}(1:14)=='## Time Closed'
                        pValue=find(header{i}==' ',7,'first');
                        try
                        obj.endDate=datenum(header{i}([pValue(4)+1:pValue(5) pValue(7)+1:end]));
                        catch
                              obj.endDate='File was not closed properly';
                        end
                    end
                end
            end
            obj.sample_us=1/obj.samplingFrequency*1e6;
            obj.saturationVoltage_uV=obj.A2DMaxValue*obj.A2DBitVolts;
            obj.recordLength_us=obj.samplesPerRecord*obj.sample_us;
            
            %detect gap in timestamp, take the last stretch of continuous
            %timestamps as the dataset.
            timeStampGap=any(find(diff(timeStamps)~=obj.recordLength_us));
            if timeStampGap
                disp('Warning, discontinuous timestamps detected. Taking the last continuous stretch as the dataset')
                gapEnd=max(find(diff(timeStamps)~=obj.recordLength_us))+1;
                obj.firstRecordStartTime_us=timeStamps(gapEnd);
                obj.lastRecordStartTime_us=timeStamps(end);
            else
                obj.firstRecordStartTime_us=timeStamps(1);
                obj.lastRecordStartTime_us=timeStamps(end);
            end
            
            obj.recordingDuration_ms=(obj.lastRecordStartTime_us-obj.firstRecordStartTime_us+obj.recordLength_us)/1000;
        end
    end
end