classdef MCRackRecording < dataRecording
    properties
        %{
        Properties of parent class: dataRecording
        recordingName=[]; %The name of the recording
        recordingDir %Full directory containing the recorded session
        dataFileNames %Name of recording data files
        startDate %Start date (time) of Recording (matlab long format)
        endDate %End date (time) of Recording (matlab long format)
        samplingFrequency %Sampling rate
        recordingDuration_ms %the total duration of the recording in [ms]
        channelNames %a cell array with the names of the channels
        channelNumbers %an array with integer channel numbers
        triggerNames %the names of trigger channels
        dspLowCutFrequency %Low-pass cutoff frequency in the Neuralynx DSP (in raw data)
        dspHighCutFrequency %High-pass cutoff frequency in the Neuralynx DSP (in raw data)
        multifileMode %if multi files were selected
        nRecordings %number of recordings
        %}
        includeDigitalDataInTriggers=1;
        includeOnlyDigitalDataInTriggers=0;
        maxNumberOfDigitalChannels=4;
    end
    properties (Constant, Hidden)
        defaultRawDataStreamName='Electro';
        defaultFilteredDataStreamName='Filtere';
        defaultAnalogDataStreamName='Analog '; %dont erase the space
        defaultTriggerStreamName='Trigger';
        defaultDigitalDataStreamName='Digital';
        fileExtension='mcd';
    end
    properties (Constant)
        defaultLocalDir='C:\Users\Tulip\Documents\Academic\Post-Doc\Experiments'; %Default directory from which search starts
        signalBits = 16; %the quantization of the sampling card
        numberOfCharFromEndToBaseName=7;
    end
    properties (Hidden)
        MCinfo
        getDataConfig
        totalChannels
        totalAnalogChannels
        streamNames
        rawDataStreamNumber
        rawDataInfo
        filteredDataStreamNumber
        analogDataStreamNumber
        analogDataInfo
        triggerStreamNumber
        triggerInfo
        digitalDataStreamNumber
        digitalDataInfo
        channelID
        number2ID
        ID2number
        ZeroADValueAnalog
        MicrovoltsPerADAnalog
        SweepStartTime
        SweepStopTime
        fileOpenStruct
        recordingDurationLocal_ms
        startDateLocal
        endDateLocal
        cumEnd
        cumStart
    end
    methods
        
        function [V_uV,T_ms]=getData(obj,channels,startTime_ms,window_ms,name)
            %Extract MCRack recording raw data from file
            %Usage: [V_uV,T_ms]=obj.getData(channels,startTime_ms,window_ms,name);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the stream (if not entered, default name is used)
            %Output: V_uV - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        T_ms - A time vector relative to recording start (t=0 at start)
            if nargin == 2
                startTime_ms=0;
                window_ms=obj.recordingDuration_ms;
                obj.getDataConfig.streamname=obj.streamNames{obj.rawDataStreamNumber};
            elseif nargin==4
                obj.getDataConfig.streamname=obj.streamNames{obj.rawDataStreamNumber};
            elseif nargin==5
                obj.getDataConfig.streamname=name;
                %this option should be revised because currently all parameters are derived from the raw data stream
            else
                error('wrong number of inputs,use: [V_uV,T_ms]=getData(obj,channels,startTime_ms,window_ms,name)');
            end
            if isempty(channels) %if no channels are entered, get all channels
                channels=obj.channelNumbers;
            end
            conversionFactor=1/1000*obj.samplingFrequency;
            startTime_ms=round(startTime_ms*conversionFactor)/conversionFactor;
            window_ms=round(window_ms*conversionFactor)/conversionFactor;
            endTime_ms=startTime_ms+window_ms; %no need to conversion factor
            recordingDuration_ms=round(obj.recordingDuration_ms*conversionFactor)/conversionFactor;
            windowSamples=round(window_ms*conversionFactor);
            nTrials=length(startTime_ms);
            V_uV=ones(numel(channels),nTrials,windowSamples,'uint16')*obj.ZeroADValue;
            
            cumStart=[-Inf obj.cumStart obj.cumEnd(end)];
            cumEnd=[0 obj.cumEnd Inf];
            obj.getDataConfig.StreamNumber=obj.rawDataStreamNumber-1;
            if obj.multifileMode %this mode currently does not support extraction from edges of the recording
                for i=1:nTrials
                    tmpStartTime=startTime_ms(i);
                    startSample=1;
                    
                    pFileStart=find(startTime_ms(i)>=cumStart,1,'last');
                    pFileEnd=find((startTime_ms(i)+window_ms)<=cumEnd,1,'first');
                    
                    for f=pFileStart:pFileEnd
                        tmpEndTime=min([cumEnd(f) endTime_ms(i)]);
                        endSample=round(startSample+(tmpEndTime-tmpStartTime)/1000*obj.samplingFrequency)-1;
                        
                        if f>1 && f<=(obj.nRecordings+1) % data in inside recording range
                            mcstreammex(obj.fileOpenStruct(f-1));
                            obj.getDataConfig.startend=[tmpStartTime;tmpEndTime]-cumStart(f);
                            data=mcstreammex(obj.getDataConfig);
                            data=reshape(data.data,obj.totalChannels,length(data.data)/obj.totalChannels);
                            V_uV(:,i,startSample:endSample)=data(obj.number2ID(channels),:);
                        else % some of the data is outside the recording range - add zeros
                            V_uV(:,i,startSample:endSample)=obj.ZeroADValue;
                        end
                        startSample=endSample+1;
                        tmpStartTime=tmpEndTime;
                    end
                end
            else
                for i=1:nTrials
                    obj.getDataConfig.startend=[startTime_ms(i);startTime_ms(i)+window_ms];
                    if startTime_ms(i)>=0 && (startTime_ms(i)+window_ms)<=recordingDuration_ms
                        data=mcstreammex(obj.getDataConfig);
                        data=reshape(data.data,obj.totalChannels,length(data.data)/obj.totalChannels);
                        V_uV(:,i,:)=data(obj.number2ID(channels),:);
                    else
                        startSample=min(0,round(startTime_ms(i)*conversionFactor));
                        endSample=min(windowSamples,round((recordingDuration_ms-startTime_ms(i))*conversionFactor)); %end sample in window (not in recroding)
                        obj.getDataConfig.startend=[max(0,startTime_ms(i));min(startTime_ms(i)+window_ms,recordingDuration_ms)];
                        data=mcstreammex(obj.getDataConfig);
                        data=reshape(data.data,obj.totalChannels,length(data.data)/obj.totalChannels);
                        V_uV(:,i,1-startSample:endSample)=data(obj.number2ID(channels),:);
                        disp('Recording at edge');
                    end
                end
            end
            if obj.convertData2Double
                V_uV = (double(V_uV) - obj.ZeroADValue) * obj.MicrovoltsPerAD;
            end
            
            if nargout==2
                T_ms=(1:windowSamples)*(1e3/obj.samplingFrequency);
            end
        end
        
        function [V_uV,T_ms]=getAnalogData(obj,channels,startTime_ms,window_ms,name)
            %Extract MCRack recording analog data
            %Usage: [V_uV,T_ms]=obj.getAnalogData(channels,startTime_ms,window_ms,name);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the stream (if not entered, default name is used)
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        T_ms - A time vector relative to recording start (t=0 at start)
            
            %obj.analogDataStreamNumber=find( cellfun(@(x) all(x(1:numel(obj.defaultAnalogDataStreamName))==obj.defaultAnalogDataStreamName),obj.streamNames) );
                %obj.analogDataInfo
            if nargin==4
                obj.getDataConfig.streamname=obj.streamNames{obj.analogDataStreamNumber};
            elseif nargin==5
                obj.getDataConfig.streamname=name;
                %this option should be revised because currently all parameters are derived from the raw data stream
            else
                error('method getAnalogData was not used correctly: wrong number of inputs');
            end
            if isempty(channels)
                channels=obj.analogChannelNumbers;
            end
            conversionFactor=1/1000*obj.samplingFrequency;
            startTime_ms=round(startTime_ms*conversionFactor)/conversionFactor;
            window_ms=round(window_ms*conversionFactor)/conversionFactor;
            endTime_ms=startTime_ms+window_ms; %no need to conversion factor
            recordingDuration_ms=round(obj.recordingDuration_ms*conversionFactor)/conversionFactor;
            windowSamples=round(window_ms*conversionFactor);
            nTrials=length(startTime_ms);
            V_uV=ones(numel(channels),nTrials,windowSamples,'uint16')*obj.ZeroADValueAnalog;
            
            cumStart=[-Inf obj.cumStart obj.cumEnd(end)];
            cumEnd=[0 obj.cumEnd Inf];
            obj.getDataConfig.StreamNumber=obj.analogDataStreamNumber-1;
            if obj.multifileMode %this mode currently does not support extraction from edges of the recording
                for i=1:nTrials
                    tmpStartTime=startTime_ms(i);
                    startSample=1;
                    
                    pFileStart=find(startTime_ms(i)>=cumStart,1,'last');
                    pFileEnd=find((startTime_ms(i)+window_ms)<=cumEnd,1,'first');
                    
                    for f=pFileStart:pFileEnd
                        tmpEndTime=min([cumEnd(f) endTime_ms(i)]);
                        endSample=round(startSample+(tmpEndTime-tmpStartTime)/1000*obj.samplingFrequency)-1;
                        
                        if f>1 && f<=(obj.nRecordings+1) % data in inside recording range
                            mcstreammex(obj.fileOpenStruct(f-1));
                            obj.getDataConfig.startend=[tmpStartTime;tmpEndTime]-cumStart(f);
                            data=mcstreammex(obj.getDataConfig);
                            data=reshape(data.data,obj.totalAnalogChannels,length(data.data)/obj.totalAnalogChannels);
                            V_uV(:,i,startSample:endSample)=data(channels,:);
                        else % some of the data is outside the recording range - add zeros
                            V_uV(:,i,startSample:endSample)=obj.ZeroADValueAnalog;
                        end
                        startSample=endSample+1;
                        tmpStartTime=tmpEndTime;
                    end
                end
            else
                for i=1:nTrials
                    obj.getDataConfig.startend=[startTime_ms(i);startTime_ms(i)+window_ms];
                    if startTime_ms(i)>=0 && (startTime_ms(i)+window_ms)<=recordingDuration_ms
                        data=mcstreammex(obj.getDataConfig);
                        data=reshape(data.data,obj.totalAnalogChannels,length(data.data)/obj.totalAnalogChannels);
                        V_uV(:,i,:)=data(channels,:);
                    else
                        startSample=min(0,round(startTime_ms(i)*conversionFactor));
                        endSample=min(windowSamples,round((recordingDuration_ms-startTime_ms(i))*conversionFactor)); %end sample in window (not in recroding)
                        obj.getDataConfig.startend=[max(0,startTime_ms(i));min(startTime_ms(i)+window_ms,recordingDuration_ms)];
                        data=mcstreammex(obj.getDataConfig);
                        data=reshape(data.data,obj.totalAnalogChannels,length(data.data)/obj.totalAnalogChannels);
                        V_uV(:,i,1-startSample:endSample)=data(obj.number2ID(channels),:);
                        disp('Recording at edge');
                    end
                end
            end
            if obj.convertData2Double
                V_uV = (double(V_uV) - obj.ZeroADValueAnalog) * obj.MicrovoltsPerADAnalog;
            end
            
            if nargout==2
                T_ms=(1:windowSamples)*(1e3/obj.samplingFrequency);
            end
        end
        
        function [D,T_ms]=getDigitalData(obj,startTime_ms,window_ms,name)
            %Extract MCRack digital data from recording
            %Usage: [D,T_ms]=getDigitalData(startTime_ms,window_ms,name)
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the stream (if not entered, default name is used)
            %Output: D - A 3D matrix [nChannels x nTrials x nSamples] with digitalData waveforms across specified channels and trials
            %        T_ms - A time vector relative to recording start (t=0 at start)
            if nargin==3
                obj.getDataConfig.streamname=obj.streamNames{obj.digitalDataStreamNumber};
            elseif nargin==4
                obj.getDataConfig.streamname=name;
                %this option should be revised because currently all parameters are derived from the raw data stream
            elseif nargin==1
                startTime_ms=0;
                window_ms=obj.recordingDuration_ms;
                obj.getDataConfig.streamname=obj.streamNames{obj.digitalDataStreamNumber};
            else
                error('method getData was not used correctly: wrong number of inputs');
            end
            conversionFactor=1/1000*obj.samplingFrequency;
            startTime_ms=round(startTime_ms*conversionFactor)/conversionFactor;
            window_ms=round(window_ms*conversionFactor)/conversionFactor;
            endTime_ms=startTime_ms+window_ms; %no need to conversion factor
            recordingDuration_ms=round(obj.recordingDuration_ms*conversionFactor)/conversionFactor;
            windowSamples=round(window_ms*conversionFactor);
            nTrials=length(startTime_ms);
            
            D=false(obj.maxNumberOfDigitalChannels,nTrials,windowSamples); %up to 4 digital bits are allowed
            
            obj.getDataConfig.StreamNumber=obj.digitalDataStreamNumber-1;
            if obj.multifileMode %this mode currently does not support extraction from edges of the recording
                for i=1:nTrials
                    pFileStart=find(startTime_ms(i)>=obj.cumStart,1,'last');
                    pFileEnd=find((startTime_ms(i)+window_ms)<=obj.cumEnd,1,'first');
                    tmpStartTime=startTime_ms(i);
                    startSample=1;
                    for f=pFileStart:pFileEnd
                        mcstreammex(obj.fileOpenStruct(f));
                        
                        tmpEndTime=min([obj.cumEnd(f) endTime_ms(i)]);
                        obj.getDataConfig.startend=[tmpStartTime;tmpEndTime]-obj.cumStart(f);
                        
                        if tmpStartTime>=0 && tmpEndTime<=recordingDuration_ms
                            data=mcstreammex(obj.getDataConfig);
                            endSample=startSample+numel(data.data)-1;
                            D(:,i,startSample:endSample)=rem(floor(data.data*pow2(0:-1:(1-obj.maxNumberOfDigitalChannels))),2)';
                        else
                            error('Requested data is outside stream limits - this is currently not supported in multi file mode');
                        end
                        startSample=endSample+1;
                        tmpStartTime=tmpEndTime;
                    end
                end
            else %single file mode
                for i=1:nTrials
                    obj.getDataConfig.startend=[startTime_ms(i);startTime_ms(i)+window_ms];
                    if startTime_ms(i)>=0 && (startTime_ms(i)+window_ms)<=recordingDuration_ms
                        data=mcstreammex(obj.getDataConfig);
                        D(:,i,:)=rem(floor(data.data*pow2(0:-1:(1-obj.maxNumberOfDigitalChannels))),2)';
                    else
                        startSample=min(0,round(startTime_ms(i)*conversionFactor));
                        endSample=min(windowSamples,round((recordingDuration_ms-startTime_ms)*conversionFactor));
                        obj.getDataConfig.startend=[max(0,startTime_ms(i));min(startTime_ms(i)+window_ms,recordingDuration_ms)];
                        data=mcstreammex(obj.getDataConfig);
                        D(:,i,1-startSample:endSample)=rem(floor(data.data*pow2(0:-1:(1-obj.maxNumberOfDigitalChannels))),2)';
                        disp('Recording at edge');
                    end
                end
            end
            if nargout==2
                T_ms=(1:windowSamples)*(1e3/obj.samplingFrequency);
            end
        end
        
        function [T_ms]=getTrigger(obj,startTime_ms,window_ms,name)
            %Extract triggers from recording
            %Usage : [T_ms]=obj.getTrigger(,startTime_ms,window_ms,name)
            %Input : startTime_ms - start time [ms].
            %        window_ms - the window duration [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        name
            %Output: T_ms - trigger times [ms] - different triggers are arranged in a cell array
            if isempty(obj.triggerStreamNumber) && isempty(obj.digitalDataStreamNumber) 
                disp('Warming!!!! No trigger or digital data entities in recording');
                T_ms={};
                return;
            end
            if nargin==4
                obj.getDataConfig.streamname=name;
            elseif nargin==3
                obj.getDataConfig.streamname=obj.streamNames(obj.triggerStreamNumber);
            elseif nargin==1
                startTime_ms=0;
                window_ms=obj.recordingDuration_ms;
                obj.getDataConfig.streamname=obj.streamNames(obj.triggerStreamNumber);
            end
            
            endTime_ms=startTime_ms+window_ms; %no need to conversion factor
            
            pFileStart=find(startTime_ms>=obj.cumStart,1,'last');
            pFileEnd=find((startTime_ms+window_ms)<=obj.cumEnd,1,'first');
            nFiles=pFileEnd-pFileStart+1;
            
            if obj.includeOnlyDigitalDataInTriggers
                nTriggers=0;
            else
                nTriggers=numel(obj.triggerStreamNumber);
            end
            nTriggersDigital=2*obj.maxNumberOfDigitalChannels;
            T=cell(nTriggers+nTriggersDigital,nFiles);
            
            tmpStartTime=startTime_ms;
            for i=pFileStart:pFileEnd
                tmpEndTime=min([obj.cumEnd(i) endTime_ms]);
                        
                obj.MCinfo=mcstreammex(obj.fileOpenStruct(i)); %a file can be read only after opening, and if not other file was opened afterwards
                obj.getDataConfig.startend=[tmpStartTime;tmpEndTime]-obj.cumStart(i);
                
                for j=1:nTriggers
                    obj.getDataConfig.StreamNumber=obj.triggerStreamNumber-1;
                    triggerData=mcstreammex(obj.getDataConfig);
                    T{j,i}=triggerData.times+obj.cumStart(i);
                end
                
                if obj.includeDigitalDataInTriggers & ~isempty(obj.digitalDataStreamNumber)
                    [D,Ttmp]=getDigitalData(obj,tmpStartTime,tmpEndTime-tmpStartTime);
                    %validChannels=any(any(D,3),2); %problematic since trigger identity changes
                    %D=D(validChannels,:,:);
                    for j=1:size(D,1)
                        T{nTriggers+2*j-1,i}=Ttmp(find(diff(squeeze(D(j,:,:)))>0))+obj.cumStart(i);
                        T{nTriggers+2*j,i}=Ttmp(find(diff(squeeze(D(j,:,:)))<0))+obj.cumStart(i);
                    end
                end
                tmpStartTime=tmpEndTime;
            end
            
            for i=1:(nTriggers+nTriggersDigital)
                T_ms{i}=cell2mat(T(i,:));
            end
            %T_ms(cellfun(@isempty,T_ms))=[];
            
        end
    end
    
    methods (Hidden)
        %class constructor
        function obj = MCRackRecording(recordingFile)
            %Usage: obj = MCRackRecording(recordingFile),
            %structure of resulting file: [nChannels x nTrials x nSamples]

            obj.datatype='uint16';

            %get data files
            if nargin==0
                recordingFile=[];
            elseif nargin>1
                disp('Object was not constructed since too many parameters were given at construction');
                return;
            end
            obj=obj.getRecordingFiles(recordingFile,obj.fileExtension);
            
            %open files for reading and getting file data
            for i=1:obj.nRecordings
                tmpFileName=[obj.recordingDir filesep obj.dataFileNames{i}];
                obj.fileOpenStruct(i).function='OpenFile';
                if ~exist(tmpFileName,'file')
                    error(['Recording file does not exist: ' tmpFileName]);
                end
                obj.fileOpenStruct(i).Filename=[obj.recordingDir filesep obj.dataFileNames{i}];
                obj.MCinfo=mcstreammex(obj.fileOpenStruct(i)); %a file can be read only after opening, and if no other file set was opened afterwards
                
                obj.SweepStartTime(i)=obj.MCinfo.SweepStartTime(1);
                obj.SweepStopTime(i)=obj.MCinfo.SweepStopTime(end);
                
                %notice that these fields are not the same as the ones in data tool (2 hour delay)
                obj.startDateLocal(i)=datenum(obj.MCinfo.RecordingDate);
                obj.endDateLocal(i)=datenum(obj.MCinfo.RecordingStopDate);
                
                obj.recordingDurationLocal_ms(i) = obj.SweepStopTime(i)-obj.SweepStartTime(i);
            end
            obj.recordingDuration_ms=sum(obj.recordingDurationLocal_ms);
            obj.cumEnd=cumsum(obj.recordingDurationLocal_ms);
            obj.cumStart=[obj.SweepStartTime(1) obj.cumEnd(1:end-1)];
            obj.startDate=obj.startDateLocal(1);
            obj.endDate=obj.endDateLocal(end);
            
            if obj.nRecordings>1
                obj.multifileMode=true;
            else
                obj.multifileMode=false;
            end
            %check that files were entered in the correct order
            if obj.multifileMode
                for i=2:obj.nRecordings
                    if obj.startDateLocal(i)<obj.startDateLocal(i-1) 
                        disp('Error: mcd files were not entered in their chronological order!!! Please reorder before running');
                        return;
                    end
                end
            end
            obj.streamNames=obj.MCinfo.StreamNames;
            
            obj.rawDataStreamNumber=find( cellfun(@(x) all(x(1:numel(obj.defaultRawDataStreamName))==obj.defaultRawDataStreamName),obj.streamNames) );
            obj.filteredDataStreamNumber=find( cellfun(@(x) all(x(1:numel(obj.defaultFilteredDataStreamName))==obj.defaultFilteredDataStreamName),obj.streamNames) );
            
            if isempty(obj.rawDataStreamNumber) & ~isempty(obj.filteredDataStreamNumber) %in case there is no raw data, but a stream of filtered data exist, treat that stream as raw
                disp('Filtered data stream is taken instead of raw data');
                obj.rawDataStreamNumber=obj.filteredDataStreamNumber;
            end
            
            if ~isempty(obj.rawDataStreamNumber)
                obj.rawDataInfo=obj.MCinfo.StreamInfo{obj.rawDataStreamNumber};
                
                %This assumes that all streams had the same sampling frequency
                obj.samplingFrequency=obj.MCinfo.MillisamplesPerSecond2(obj.rawDataStreamNumber)/1000;   %in Hz
                
                obj.totalChannels=numel(find(obj.MCinfo.ChannelID(:,obj.rawDataStreamNumber)~=-1));
                
                obj=obj.loadChLayout;
                if isempty(obj.chLayoutNumbers)
                    if obj.totalChannels<=32
                        obj.layoutName='layout_300_6x6_FlexMEA';
                        obj.electrodePitch=300;
                    elseif obj.totalChannels<=60
                        obj.layoutName='layout_200_8x8.mat';
                        obj.electrodePitch=200;
                    elseif obj.totalChannels<=120
                        obj.layoutName='layout_100_12x12.mat';
                        obj.electrodePitch=100;
                    elseif obj.totalChannels<=252
                        obj.layoutName='layout_100_16x16.mat';
                        obj.electrodePitch=100;
                    end
                    load(obj.layoutName);
                    obj.chLayoutNumbers=En;
                    obj.chLayoutNames=Ena;
                    obj.chLayoutPositions=Enp;
                end
                
                obj.channelNames=obj.MCinfo.ChannelNames(:,obj.rawDataStreamNumber);
                obj.channelNumbers=obj.MCinfo.ChannelID(:,obj.rawDataStreamNumber)+1;
                
                %check that channels are in consecutive order
                if any(sort(obj.channelNumbers)~=obj.channelNumbers)
                    error('Channels do not appear in consecutive order!!!');
                end
                
                obj.channelNumbers=1:numel(obj.channelNumbers);
                
                %remove from layout channels not in the channel list
                unrecordedChannels=setdiff(obj.chLayoutNumbers(~isnan(obj.chLayoutNumbers)),obj.channelNumbers);
                for i=1:numel(unrecordedChannels)
                    p=find(obj.chLayoutNumbers==unrecordedChannels(i));
                    obj.chLayoutNumbers(p)=NaN;
                    obj.chLayoutNumbers(obj.chLayoutNumbers>unrecordedChannels(i))=obj.chLayoutNumbers(obj.chLayoutNumbers>unrecordedChannels(i))-1;
                    obj.channelNumbers(obj.channelNumbers>unrecordedChannels(i))=obj.channelNumbers(obj.channelNumbers>unrecordedChannels(i))-1;
                    if ~isempty(obj.chLayoutPositions)
                        obj.chLayoutPositions(:,unrecordedChannels(i))=0;
                    end
                    obj.chLayoutNames{p}=[];
                end
                
                if max(obj.channelNumbers)~=numel(obj.channelNumbers) || min(obj.channelNumbers)~=1
                    error('Channel numbers in MCRackRecording object do not start from one or non-continuous');
                end
                
                obj.channelID=1:obj.totalChannels;
                obj.number2ID(obj.channelNumbers)=obj.channelID;
                obj.ID2number(obj.channelID)=obj.channelNumbers;
                
                obj.ZeroADValue = obj.MCinfo.ZeroADValue(obj.rawDataStreamNumber);
                if obj.MCinfo.UnitSign(obj.rawDataStreamNumber)=='V'
                    obj.MicrovoltsPerAD = obj.MCinfo.UnitsPerAD(obj.rawDataStreamNumber)*1000000;
                else
                    error('Units in input data are wrong!!!')
                end
                
            end
            
            obj.analogDataStreamNumber=find( cellfun(@(x) all(x(1:numel(obj.defaultAnalogDataStreamName))==obj.defaultAnalogDataStreamName),obj.streamNames) );
            if ~isempty(obj.analogDataStreamNumber)
                obj.analogDataInfo=obj.MCinfo.StreamInfo{obj.analogDataStreamNumber};
                obj.analogChannelNumbers=obj.MCinfo.ChannelID(obj.MCinfo.ChannelID(:,obj.analogDataStreamNumber)~=-1,obj.analogDataStreamNumber);
                obj.totalAnalogChannels=numel(find(obj.MCinfo.ChannelID(:,obj.analogDataStreamNumber)~=-1));
                
                obj.ZeroADValueAnalog = obj.MCinfo.ZeroADValue(obj.analogDataStreamNumber);
                if obj.MCinfo.UnitSign(obj.analogDataStreamNumber)=='V'
                    obj.MicrovoltsPerADAnalog = obj.MCinfo.UnitsPerAD(obj.analogDataStreamNumber)*1000000;
                else
                    error('Units in input data are wrong!!!')
                end
                if isempty(obj.samplingFrequency) %if there is no raw data extract sampling frequency from analog data
                    obj.samplingFrequency=obj.MCinfo.MillisamplesPerSecond2(obj.analogDataStreamNumber)/1000;   %in Hz
                end
            end
            
            obj.triggerStreamNumber=find( cellfun(@(x) all(x(1:numel(obj.defaultTriggerStreamName))==obj.defaultTriggerStreamName),obj.streamNames) );
            if ~isempty(obj.triggerStreamNumber)
                obj.triggerInfo=[obj.MCinfo.StreamInfo{obj.triggerStreamNumber}];
            end
            
            obj.digitalDataStreamNumber=find( cellfun(@(x) all(x(1:numel(obj.defaultDigitalDataStreamName))==obj.defaultDigitalDataStreamName),obj.streamNames) );
            if ~isempty(obj.digitalDataStreamNumber)
                obj.digitalDataInfo=[obj.MCinfo.StreamInfo{obj.digitalDataStreamNumber}];
                if isempty(obj.samplingFrequency) %if there is no raw data and now analog data extract sampling frequency from digital data
                    obj.samplingFrequency=obj.MCinfo.MillisamplesPerSecond2(obj.digitalDataStreamNumber)/1000;   %in Hz
                end
            end
            
            %determine channel layout (using NeuroShare) - It seems that neuroshare always gives back the standard MEA configuration
            %{
            [ns] = ns_SetLibrary(which('nsMCDlibrary64.dll'));
            [~, hfile] = ns_OpenFile([obj.recordingDir filesep obj.dataFileNames{1}]);
            [~, fileInfo] = ns_GetFileInfo(hfile);
            [~, entityInfo] = ns_GetEntityInfo(hfile, 1:fileInfo.EntityCount);
            entityName = cellfun(@(x) x(1:4),{entityInfo.EntityLabel},'UniformOutput',0);
            rawElectrodeEnteties = find(strcmp(entityName,'elec'));
            [~, rawElectrodeInfo] = ns_GetAnalogInfo(hfile,rawElectrodeEnteties);
            xCoord=cell2mat({rawElectrodeInfo.LocationX});
            yCoord=cell2mat({rawElectrodeInfo.LocationY});
            %}

            %definitions for MCStreamMex function
            obj.getDataConfig.function='GetFromTo';
            obj.getDataConfig.timesonly=0;
            obj.getDataConfig.originorder='on';
            obj.getDataConfig.warning='on';
            obj.getDataConfig.ChunkNumber=-1;
        end
    end
end

%Making a new layout
%{
for i=1:252
    tmp=obj.channelNames{i};
    x(i)=tmp(1)-0;
    y(i)=str2num(tmp(2:end));
end
x=x-min(x)+1;
x(x>10)=x(x>10)-1;
x(x==17)=16;

Ena=cell(16,16);
En=nan(16,16);
for i=1:252
    Ena(x(i),y(i))=obj.channelNames(i);
    En(x(i),y(i))=obj.channelNumbers(i);
end

Ena=Ena';
En=En';
%}