classdef filteredTimeSeriesMatlab < dataRecording
    properties
        pointerToMatFile
        fileInfo
        channel2Place
        triggerNames=[];
    end
    properties (Constant, Hidden)
        defaultLocalDir='C:\Users\Tulip\Documents\Academic\Post-Doc\Experiments'; %Default directory from which search starts
    end
    properties (Hidden)

    end
    methods
        function obj=reloadMetaData(obj)
            %Reloads the meta data from directory
            %Usage: obj.reloadMetaData();
            %Input: []
            %obj.fileInfo = h5info([obj.recordingDir '\' obj.dataFileNames]);
            obj.pointerToMatFile=matfile([obj.recordingDir '\' obj.dataFileNames]);
            
            obj.startDate=obj.pointerToMatFile.startDate; %Start date (time) of Recording (matlab long format)
            obj.endDate=obj.pointerToMatFile.endDate; %End date (time) of Recording (matlab long format)
            obj.samplingFrequency=obj.pointerToMatFile.filteredSamplingFrequency; %Sampling rate
            
            obj.recordingDuration_ms=obj.pointerToMatFile.recordingDuration_ms; %the total duration of the recording in [ms]
            obj.channelNames=obj.pointerToMatFile.channelNames; %a cell array with the names of the channels
            obj.channelNumbers=obj.pointerToMatFile.channelNumbers; %an array with integer channel numbers
            obj.dspLowCutFrequency=obj.pointerToMatFile.cutoffLow; %Low-pass cutoff frequency in the Neuralynx DSP (in raw data)
            obj.dspHighCutFrequency=obj.pointerToMatFile.cutoffHigh; %High-pass cutoff frequency in the Neuralynx DSP (in raw data)
            
            obj.channel2Place(obj.channelNumbers)=obj.channelNumbers;
            TriggerVars=whos(obj.pointerToMatFile,'Trigger*');
            for i=1:numel(TriggerVars)
                obj.triggerNames{i}=TriggerVars(i).name;
            end
        end
        function [V_uV,t_ms]=getData(obj,channels,startTime_ms,window_ms)
            %Extract filtered matlab data from file to memory
            %Usage: [V_uV,t_ms]=obj.getData(channels,startTime_ms,window_ms);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials  
            %        t_ms - A time vector relative to recording start (t=0 at start)
            
            nSamples=round(window_ms/1e3*obj.samplingFrequency);
            nTrials=numel(startTime_ms);
            nCh=numel(channels);
            
            startSample=floor(startTime_ms/1e3*obj.samplingFrequency);
            totalSamplesInRecording=size(obj.pointerToMatFile,'Data',2);
            
            V_uV=NaN(nCh,nTrials,nSamples); %initialize waveform matrix
            if all(startTime_ms>=0) && all((startTime_ms+window_ms)<=(obj.recordingDuration_ms)) %fast extaction does not work on recording edges
                for i=1:nTrials
                    V_uV(:,i,:)=obj.pointerToMatFile.Data(obj.channel2Place(channels),(startSample(i)+1):(startSample(i)+nSamples));
                end
            else %this is for the special case that the segment cuts recording end or beginning
                
                pOutSizeOfRecording=startTime_ms<=(-window_ms) | startTime_ms>=obj.recordingDuration_ms;
                pOnRecordingEdge=((startTime_ms<0) | (startTime_ms+window_ms)>obj.recordingDuration_ms) & ~pOutSizeOfRecording;
                pInsideRecording=~(pOutSizeOfRecording | pOnRecordingEdge);
                
                
                %extract voltage regions that do not require special treatment
                if any(pInsideRecording)
                    [V_uV(:,pInsideRecording,:)]=obj.getData(obj.channel2Place(channels),startTime_ms(pInsideRecording),window_ms);
                end
                
                %Go over each trial and channel and extract
                for i=find(pOnRecordingEdge)
                    startInd=max(startSample(i),1);
                    endIdx=min(startInd+nSamples-1,totalSamplesInRecording);
                    V_uV(:,i,(startInd-startSample(i)+1):(endIdx-(startInd+nSamples-1)+nSamples))=obj.pointerToMatFile.Data(obj.channel2Place(channels),startInd:endIdx);
                end
            end
            if nargout==2
                t_ms=(1:nSamples)*(1e3/obj.samplingFrequency);
            end
        end
        function [T_ms]=getTrigger(obj,startTime_ms,window_ms,name)
            %Extract triggers from file Neuralynx recording
            %Usage : [T_ms]=obj.getTrigger(startTime_ms,endTime_ms,TTLbits)
            %Input : startTime_ms - start time [ms].
            %        window_ms - the window duration [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        TTLbits - which bit to extract for time stamps (out of 8,default = first bit, 1)
            %Output: T_ms - trigger times [ms]
            if isempty(obj.triggerNames)
                disp('No triggers in recording');
                return;
            end
            if nargin==1
                startTime_ms=[];
                window_ms=[];
                name=[];
            elseif nargin==3
                name=[];
            elseif nargin>4 || nargin==2
                error('Incorrect number of inputs: [T_ms]=obj.getTrigger(startTime_ms,endTime_ms,triggerName)');
            end

            if isempty(name)
                name=obj.triggerNames;
            end
            for i=1:numel(name)
                p=strcmp(name{i},obj.triggerNames);
                T_ms{i}=obj.pointerToMatFile.(obj.triggerNames{p});
            end
        end
    end
    
    methods (Hidden)
        %class constructor
        function obj = filteredTimeSeriesMatlab(recordingFile)
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
                [obj.dataFileNames,obj.recordingDir]= uigetfile('*.mat','Choose the MCD file',obj.defaultLocalDir,'MultiSelect','on');
                if obj.dataFileNames==0 %no folder chosen
                    error('Object was not constructed since no folder was choosen');
                end
            else
                error('Object was not constructed since too many parameters were given at construction');
            end
            
            obj=reloadMetaData(obj);
        end
    end
end