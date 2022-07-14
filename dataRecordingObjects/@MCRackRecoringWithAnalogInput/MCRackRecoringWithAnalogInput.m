classdef MCRackRecoringWithAnalogInput < dataRecording
    properties
        dataFile
        streamNames
        rawDataStreamNumber
        rawDataInfo
        triggerStreamNumbers
        triggerInfo
        analogStreamNumbers
        analogInfo
        MEAType
        analogChannelNames
        electrodeChannelNames
        analogChannelID
        electrodeChannelID
        electrodeChannelNumbers
        analogChannelNumbers
        channelID
        number2ID
        ID2number
        %{
        Properties of parent class: dataRecording
        recordingDir %Full directory containing the recorded session
        startDate %Start date (time) of Recording (matlab long format)
        endDate %End date (time) of Recording (matlab long format)
        samplingFrequency %Sampling rate
        dataFileNames %Name of recording data files
        channelNames %a cell array with the names of the channels
        channelNumbers %an array with integer channel numbers
        recordingDuration_ms %the total duration of the recording
        %}
    end
    properties (Constant, Hidden)
        defaultRawDataStreamName='Electrode Raw Data';
        defaultTriggerStreamName='Trigger';
        defaultAnalogStreamName='Analog'; %'Analog Raw Data'
        defaultLocalDir='C:\Users\Tulip\Documents\Academic\Post-Doc\Experiments'; %Default directory from which search starts
    end
    properties (Hidden)
        
    end
    methods
        function [V_uV,t_ms]=getData(obj,channels,startTime_ms,window_ms,name)
            %Extract Neuralynx recording data from file to memory
            %Usage: [V_uV,t_ms]=obj.getData(channels,startTime_ms,window_ms);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        t_ms - A time vector relative to recording start (t=0 at start)
            if nargin==4
                streamName=obj.streamNames{obj.rawDataStreamNumber};
                streamNameAnalog=obj.streamNames{obj.analogStreamNumbers};
            elseif nargin==5
                streamName=name;
                streamNameAnalog=[];
            else 
                error('method getData was not used correctly: wrong number of inputs');
            end
            startTime_ms=round(startTime_ms*10)/10;
            window_ms=round(window_ms*10)/10;
            windowSamples=round(window_ms/1000*obj.samplingFrequency);
            nTrials=length(startTime_ms);
            V_uV=nan(numel(channels),nTrials,windowSamples);
            
            for i=1:nTrials
                startSample=round(startTime_ms(i)/1000*obj.samplingFrequency);
                if startSample>=0
                    electrodeData = nextdata(obj.dataFile,'streamname',streamName,'startend',[startTime_ms;startTime_ms+window_ms]);
                    analogData = nextdata(obj.dataFile,'streamname',streamNameAnalog,'originorder','on','startend',[startTime_ms;startTime_ms+window_ms]);
                    data=[ad2muvolt(obj.dataFile,electrodeData.data,streamName);...
                        ad2muvolt(obj.dataFile,analogData.data,streamNameAnalog)];
                    V_uV(:,i,:)=data(obj.number2ID(channels),:);
                else
                    electrodeData = nextdata(obj.dataFile,'streamname',streamName,'startend',[0;startTime_ms+window_ms]);
                    analogData = nextdata(obj.dataFile,'streamname',streamNameAnalog,'originorder','on','startend',[0;startTime_ms+window_ms]);
                    data=[ad2muvolt(obj.dataFile,electrodeData.data,streamName);...
                        ad2muvolt(obj.dataFile,analogData.data,streamNameAnalog)];
                    V_uV(:,i,1-startSample:end)=data(obj.number2ID(channels),:);
                end
            end
            if nargout==2
                t_ms=(1:windowSamples)*(1e3/obj.samplingFrequency);
            end
        end
        function [T_ms]=getTrigger(obj,startTime_ms,window_ms,name)
            %Extract triggers from file Neuralynx recording
            %Usage : [T_ms]=obj.getTrigger(startTime_ms,endTime_ms,TTLbits)
            %Input : startTime_ms - start time [ms].
            %        window_ms - the window duration [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        TTLbits - which bit to extract for time stamps (out of 8,default = first bit, 1)
            %Output: T_ms - trigger times [ms] - different triggers are arranged in a cell array
            if isempty(obj.triggerStreamNumbers)
                disp('Warming!!!! No trigger entities in recording');
                T_ms=NaN;
                return;
            end
            if nargin==4
                streamName=name;
            elseif nargin==3
                streamName=obj.streamNames(obj.triggerStreamNumbers);
            elseif nargin==1
                startTime_ms=0;
                window_ms=obj.recordingDuration_ms;
                streamName=obj.streamNames(obj.triggerStreamNumbers);
            end
            for i=1:numel(obj.triggerStreamNumbers)
                c = nextdata(obj.dataFile,'streamname',streamName{i},'startend',[startTime_ms;startTime_ms+window_ms]);
                T_ms{i}=c.times;
            end
        end
    end
    
    methods (Hidden)
        %class constructor
        function obj = MCRackRecoringWithAnalogInput(recordingFile)
            if nargin==1 %if directory with data was not entered open get directory GUI
                [pathstr, name, ext] = fileparts(recordingFile);
                if isempty(pathstr) %in case the file is in the current directory
                    obj.recordingDir=cd;
                else
                    obj.recordingDir=pathstr;
                end
                obj.dataFileNames=[name ext];
                if ~isdir(obj.recordingDir) || ~exist([obj.recordingDir '\' obj.dataFileNames],'file')
                    error('Object was not constructed since no valid folder was choosen');
                end
            elseif nargin==0
                [obj.dataFileNames,obj.recordingDir]= uigetfile('*.mcd','Choose the MCD file',obj.defaultLocalDir,'MultiSelect','on');
                if obj.dataFileNames==0 %no folder chosen
                    error('Object was not constructed since no folder was choosen');
                end
            else
                error('Object was not constructed since too many parameters were given at construction');
            end
            obj.dataFile=datastrm([obj.recordingDir '\' obj.dataFileNames]);
            obj.streamNames=getfield(obj.dataFile,'StreamNames');
            obj.rawDataStreamNumber=find(strcmp(obj.streamNames,obj.defaultRawDataStreamName));
            obj.rawDataInfo=getfield(obj.dataFile,'StreamInfo',obj.rawDataStreamNumber);
            
            obj.triggerStreamNumbers=find( cellfun(@(x) all(x(1:numel(obj.defaultTriggerStreamName))==obj.defaultTriggerStreamName),obj.streamNames) );
            obj.triggerInfo=getfield(obj.dataFile,'StreamInfo',obj.triggerStreamNumbers);

            obj.analogStreamNumbers=find( cellfun(@(x) all(x(1:numel(obj.defaultAnalogStreamName))==obj.defaultAnalogStreamName),obj.streamNames) );
            obj.analogInfo=getfield(obj.dataFile,'StreamInfo',obj.analogStreamNumbers);
            
            obj.samplingFrequency=getfield(obj.dataFile,'MillisamplesPerSecond2',obj.rawDataStreamNumber)/1000;   %in Hz
            obj.startDate=getfield(obj.dataFile,'recordingdate');
            obj.endDate=getfield(obj.dataFile,'recordingStopDate');
            
            obj.recordingDuration_ms = getfield(obj.dataFile,'sweepStopTime')-getfield(obj.dataFile,'sweepStartTime');
            obj.MEAType=getfield(obj.dataFile,'meatype');
            
            %ChannelNames2           cell  (2  1)
            %HardwareChannelNames2   cell  (2  1)
            %ChannelID2              cell  (2  1)
            %HardwareChannelID2
            
            obj.electrodeChannelNames=getfield(obj.dataFile,'ChannelNames2',obj.rawDataStreamNumber);
            obj.electrodeChannelNames=obj.electrodeChannelNames{1};
            obj.electrodeChannelID=getfield(obj.dataFile,'ChannelID2',obj.rawDataStreamNumber);
            obj.electrodeChannelID=obj.electrodeChannelID{1};
            obj.electrodeChannelNumbers=str2double(obj.electrodeChannelNames);
            
            obj.analogChannelNames=getfield(obj.dataFile,'ChannelNames2',obj.analogStreamNumbers);
            obj.analogChannelNames=obj.analogChannelNames{1};
            %obj.analogChannelID=getfield(obj.dataFile,'ChannelID2',obj.analogStreamNumbers);
            %obj.analogChannelID=obj.analogChannelID{1};
            totalCh=getfield(obj.dataFile,'TotalChannels');
            obj.analogChannelID=((totalCh+1):(totalCh+numel(obj.analogChannelNames)))'; %in this version analog channels are added to regular channels and placed right after them
            obj.analogChannelNumbers=(1:numel(obj.analogChannelID))';
            
            obj.channelNumbers=[obj.electrodeChannelNumbers;obj.analogChannelNumbers];
            obj.channelNames=[obj.electrodeChannelNames;obj.analogChannelNames];
            obj.channelID=[obj.electrodeChannelID;obj.analogChannelID];
            
            obj.number2ID(obj.channelNumbers)=obj.channelID;
            obj.ID2number(obj.channelID)=obj.channelNumbers;
        end
    end
end