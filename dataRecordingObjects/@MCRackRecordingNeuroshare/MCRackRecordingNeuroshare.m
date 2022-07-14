classdef MCRackRecordingNeuroshare < dataRecording & handle
    properties
        %{
        Properties of parent class: dataRecording
        recordingDir %Full directory containing the recorded session
        startDate %Start date (time) of Recording (matlab long format)
        endDate %End date (time) of Recording (matlab long format)
        samplingFrequency %Sampling rate
        dataFileNames %Name of recording data files
        channelNames %a cell array with the names of the channels
        channelNumbers %an array with integer channel numbers
        triggerNames %the names of trigger channels
        recordingDuration_ms %the total duration of the recording
        %}
        includeDigitalDataInTriggers=4;
    end
    properties (Constant)
        defaultLocalDir='C:\Users\Tulip\Documents\Academic\Post-Doc\Experiments'; %Default directory from which search starts
        signalBits = 16; %the quantization of the sampling card
        numberOfCharFromEndToBaseName=7;
    end
    properties (Constant, Hidden)
        defaultRawDataName='elec0001';
        defaultTriggerName='trig0001';
        defaultDigitalDataName='digi0001';
        fileExtension='mcd';
        maxNumberOfDigitalChannels=4;
    end
    properties (Hidden)
        hfile
        fileInfo={};
        entityInfo={};
        entityLabels %the full label name of each entity
        channelLabel
        digitalDataChannels
        entityTypes
        entityName
        entityCount
        rawElectrodeEnteties %the entity places with raw electrode data
        rawElectrodeInfo %the entity infos with raw electrode data
        trigggerEnteties %the entity places with trigger data
        triggerInfo
        digitalDataInfo
        digitalDataEnteties
        rawChNumber2EntityNumber %translation between raw channel numbers and entities
        recordingDurationLocal_ms
        startDateLocal
        endDateLocal
        cumEnd
        cumStart
        samplesPerMS
        dynamicRange
        totalChannels
    end
    methods
        function [V_uV,t_ms]=getData(obj,channels,startTime_ms,window_ms,name)
            %Extract MCRack recording data from file to memory
            %Usage: [V_uV,T_ms]=obj.getData(channels,startTime_ms,window_ms,name);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the stream (if not entered, default name is used)
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        T_ms - A time vector relative to recording start (t=0 at start)
            conversionFactor=1/1000*obj.samplingFrequency;
            startTime_ms=round(startTime_ms*conversionFactor)/conversionFactor;
            window_ms=round(window_ms*conversionFactor)/conversionFactor;
            endTime_ms=startTime_ms+window_ms; %no need to conversion factor
            recordingDuration_ms=round(obj.recordingDuration_ms*conversionFactor)/conversionFactor;
            windowSamples=round(window_ms*conversionFactor);
            nTrials=length(startTime_ms);
            
            if isempty(channels) %if no channels are entered, get all channels
                channels=obj.channelNumbers;
            end
            V_uV=zeros(numel(channels),nTrials,windowSamples);
            
            cumStart=[-Inf obj.cumStart obj.cumEnd(end)];
            cumEnd=[0 obj.cumEnd Inf];
            if obj.multifileMode %this mode currently does not support extraction from edges of the recording
                    
                lastF=0;
                    
                for i=1:nTrials
                    tmpStartTime=startTime_ms(i);
                    startSample=1; %for sample number in V_uV
                    
                    pFileStart=find(startTime_ms(i)>=cumStart,1,'last');
                    pFileEnd=find((startTime_ms(i)+window_ms)<=cumEnd,1,'first');
                
                    for f=pFileStart:pFileEnd
                        tmpEndTime=min([cumEnd(f) endTime_ms(i)]);
                        
                        startSampleLocal=tmpStartTime*obj.samplesPerMS; %for sample number in neuro-share - like tmpStartTime
                        windowSamplesLocal=tmpEndTime*obj.samplesPerMS-startSampleLocal; %for sample number in neuro-share 
                        endSample=startSample+windowSamplesLocal-1;
                        
                        if f>1 && f<=(obj.nRecordings+1) % data in inside recording range
                          if f~=lastF  
                            %[~, obj.hfile(f-1)] = ns_OpenFile([obj.recordingDir filesep obj.dataFileNames{f-1}]); %temporary , check if works
                            %ns_CloseFile([obj.recordingDir filesep obj.dataFileNames{f-1}]); %temporary , check if works
                          end
                            [ns,count,data]=ns_GetAnalogData(obj.hfile(f-1),obj.rawChNumber2EntityNumber(channels),startSampleLocal-cumStart(f)*obj.samplesPerMS+1,windowSamplesLocal);
                            V_uV(:,i,startSample:endSample)=data'*1e6;
                        
                        else
                            V_uV(:,i,startSample:endSample)=obj.ZeroADValue;
                        end
                        startSample=endSample+1;
                        tmpStartTime=tmpEndTime;
                    end
                    lastF=f;
                end
            else
                for i=1:nTrials
                    if startTime_ms(i)>=0 && (startTime_ms(i)+window_ms)<=obj.recordingDuration_ms
                        startSampleLocal=round(startTime_ms(i)/1000*obj.samplingFrequency);
                        [ns,count,data]=ns_GetAnalogData(obj.hfile,obj.rawChNumber2EntityNumber(channels),startSampleLocal+1,windowSamples);
                        V_uV(:,i,:)=data'*1e6;
                    else
                        startSample=min(0,startTime_ms(i)*obj.samplesPerMS); %for sample number in V_uV
                        endSample=min(windowSamples,(recordingDuration_ms-startTime_ms)*obj.samplesPerMS); %for sample number in V_uV
                        
                        startSampleLocal=max(0,startTime_ms(i))*obj.samplesPerMS; %for sample number in neuro-share
                        endSampleLocal=min(startTime_ms(i)+window_ms,recordingDuration_ms)*obj.samplesPerMS; %for sample number in neuro-share
                        
                        [ns,count,data]=ns_GetAnalogData(obj.hfile,obj.rawChNumber2EntityNumber(channels),1,endSampleLocal-startSampleLocal+1);
                        V_uV(:,i,1-startSample:endSample)=data'*1e6;
                        disp('Recording at edge');
                    end
                end
            end
            
            if ~obj.convertData2Double
                V_uV = uint16((V_uV./obj.MicrovoltsPerAD+obj.ZeroADValue));
            end
            if nargout==2
                t_ms=(1:windowSamples)*(1e3/obj.samplingFrequency);
            end
        end
        
        function [D,T_ms]=getDigitalData(obj,startTime_ms,window_ms,name)
            %Extract MCRack digital data from file to memory
            %Usage: [D,T_ms]=getDigitalData(startTime_ms,window_ms,name)
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the stream (if not entered, default name is used)
            %Output: D - A 3D matrix [nChannels x nTrials x nSamples] with digitalData waveforms across specified channels and trials
            %        T_ms - A time vector relative to recording start (t=0 at start)
            if nargin==3
                pDigitalEntity=1:numel(obj.digitalDataEnteties); %take all triggers in recording
            elseif nargin==4
                error('Entering digial data name is  currently not supported for neuroshare data');
                %this option should be revised because currently all parameters are derived from the raw data stream
            elseif nargin==2
                error('Not enough input variables: window is required when entering time frame');
            elseif nargin==1
                startTime_ms=0;
                window_ms=obj.recordingDuration_ms;
                pDigitalEntity=1:numel(obj.digitalDataEnteties); %take all triggers in recording
            else
                error('method getDigitalData was not used correctly: wrong number of inputs');
            end
            conversionFactor=1/1000*obj.samplingFrequency;
            startTime_ms=round(startTime_ms*conversionFactor)/conversionFactor;
            window_ms=round(window_ms*conversionFactor)/conversionFactor;
            endTime_ms=startTime_ms+window_ms; %no need to conversion factor
            recordingDuration_ms=round(obj.recordingDuration_ms*conversionFactor)/conversionFactor;
            windowSamples=round(window_ms*conversionFactor);
            nTrials=length(startTime_ms);
            
            D=false(obj.maxNumberOfDigitalChannels,nTrials,windowSamples); %up to 4 digital bits are allowed
            
            if obj.multifileMode %this mode currently does not support extraction from edges of the recording
                for i=1:nTrials
                    pFileStart=find(startTime_ms(i)>=obj.cumStart,1,'last');
                    pFileEnd=find((startTime_ms(i)+window_ms)<=obj.cumEnd,1,'first');
                    tmpStartTime=startTime_ms(i);
                    startSample=1;
                    for f=pFileStart:pFileEnd                     
                        tmpEndTime=min([obj.cumEnd(f) endTime_ms(i)]);
                        startSampleLocal=(tmpStartTime-obj.cumStart(f))*obj.samplesPerMS; %for sample number in neuro-share
                        windowSamplesLocal=((tmpEndTime-obj.cumStart(f))*obj.samplesPerMS-startSampleLocal); %for sample number in neuro-share
                        endSample=startSample+windowSamplesLocal-1;
                        %[~, obj.hfile(f)] = ns_OpenFile([obj.recordingDir filesep obj.dataFileNames{f}]); %temporary , check if works
                        
                        if tmpStartTime>=0 && tmpEndTime<=obj.recordingDuration_ms
                            [ns,count,data]=ns_GetAnalogData(obj.hfile(f),obj.digitalDataEnteties(pDigitalEntity),startSampleLocal+1,windowSamplesLocal);
                            if ~isempty(data)
                                endSample=startSample+numel(data)-1;
                                D(:,i,startSample:endSample)=rem(floor(data*pow2(0:-1:(1-obj.maxNumberOfDigitalChannels))),2)';
                            else
                                D(:,i,startSample:endSample)=false;
                            end
                        else
                            error('Requested data is outside stream limits - this is currently not supported in multi file mode');
                        end
                        startSample=endSample+1;
                        tmpStartTime=tmpEndTime;
                    end
                end
            else %single file mode
                for i=1:nTrials
                    if startTime_ms(i)>=0 && (startTime_ms(i)+window_ms)<=obj.recordingDuration_ms
                        startSampleLocal=round(startTime_ms(i)/1000*obj.samplingFrequency);
                        [ns,count,data]=ns_GetAnalogData(obj.hfile,obj.digitalDataEnteties(pDigitalEntity),startSampleLocal+1,windowSamples);
                        if ~isempty(data)
                            D(:,i,:)=rem(floor(data*pow2(0:-1:(1-obj.maxNumberOfDigitalChannels))),2)';
                        end
                    else
                        startSample=min(0,startTime_ms(i)*obj.samplesPerMS); %for sample number in V_uV
                        endSample=min(windowSamples,(recordingDuration_ms-startTime_ms)*obj.samplesPerMS); %for sample number in V_uV
                        
                        startSampleLocal=max(0,startTime_ms(i))*obj.samplesPerMS; %for sample number in neuro-share
                        endSampleLocal=min(startTime_ms(i)+window_ms,recordingDuration_ms)*obj.samplesPerMS; %for sample number in neuro-share
                       
                        [ns,count,data]=ns_GetAnalogData(obj.hfile,obj.digitalDataEnteties(pDigitalEntity),startSampleLocal+1,endSampleLocal-startSampleLocal);
                        if ~isempty(data)
                            D(:,i,1-startSample:endSample)=rem(floor(data*pow2(0:-1:(1-obj.maxNumberOfDigitalChannels))),2)';
                        else
                            disp('No digital data found');
                        end
                        disp('Recording at edge');
                    end
                end
            end
            if nargout==2
                T_ms=(1:windowSamples)*(1e3/obj.samplingFrequency);
            end
        end
        
        function [T_ms]=getTrigger(obj,startTime_ms,window_ms,name)
            %Extract triggers from file Neuralynx recording
            %Usage : [T_ms]=obj.getTrigger(startTime_ms,window_ms,name)
            %Input : startTime_ms - start time [ms].
            %        window_ms - the window duration [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        name
            %Output: T_ms - trigger times [ms] - different triggers are arranged in a cell array
            if isempty(obj.digitalDataInfo) && isempty(obj.triggerInfo)
                disp('Warming!!!! No trigger or digital data entities in recording');
                T_ms={};
                return;
            end
            if nargin==4
                error('getting triggers by name is currently not supported for neuroshare data type');
            elseif nargin==3
                pTrigEntity=1:numel(obj.trigggerEnteties); %take all triggers in recording
            elseif nargin==1
                startTime_ms=0;
                window_ms=obj.recordingDuration_ms;
                pTrigEntity=1:numel(obj.trigggerEnteties); %take all triggers in recording
            end
            
            endTime_ms=startTime_ms+window_ms; %no need to conversion factor
            
            pFileStart=find(startTime_ms>=obj.cumStart,1,'last');
            pFileEnd=find((startTime_ms+window_ms)<=obj.cumEnd,1,'first');
            nFiles=pFileEnd-pFileStart+1;
            
            nTriggers=numel(pTrigEntity);
            nTriggersDigital=2*obj.maxNumberOfDigitalChannels;
            T=cell(nTriggers+nTriggersDigital,nFiles);
            
            tmpStartTime=startTime_ms;
            for i=pFileStart:pFileEnd
                tmpEndTime=min([obj.cumEnd(i) endTime_ms]);
                
                startEndLocal=[tmpStartTime;tmpEndTime]-obj.cumStart(i);
                for j=1:nTriggers
                    [ns,triggerData] = ns_GetEventData(obj.hfile(i),obj.trigggerEnteties(pTrigEntity(j)),1:obj.entityCount{i}{obj.trigggerEnteties(pTrigEntity(j))});
                    triggerData=triggerData*1000; %change units to ms
                    triggerData=(triggerData(triggerData>=startEndLocal(1) & triggerData<startEndLocal(2)))';
                    T{j,i}=triggerData+obj.cumStart(i);
                end
                
                if obj.includeDigitalDataInTriggers & ~isempty(obj.digitalDataInfo)
                    [D,Ttmp]=getDigitalData(obj,tmpStartTime,tmpEndTime-tmpStartTime);
                    %validChannels=any(any(D,3),2);%problematic since trigger identity changes
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
        
        function obj=closeOpenFiles(obj)
            for i=1:obj.nRecordings
                [ns] = ns_CloseFile(obj.hfile(i));
            end
        end
    end
    
    methods (Hidden)
        %class constructor
        function obj = MCRackRecordingNeuroshare(recordingFile)
            %get data files
            if nargin==0
                recordingFile=[];
            elseif nargin>1
                disp('Object was not constructed since too many parameters were given at construction');
                return;
            end
            obj=obj.getRecordingFiles(recordingFile,obj.fileExtension);
            if ispc
                [ns] = ns_SetLibrary(which('nsMCDlibrary64.dll')); %finds and opens the dll library
            elseif isunix
                [ns] = ns_SetLibrary(which('nsMCDLibrary.so')); %finds and opens the dll library
            end
            %[ns,info] = ns_GetLibraryInfo();
            
            %open files for reading and getting file data
            fprintf('Opening files for reading (/%d): ',obj.nRecordings);
            for i=1:obj.nRecordings
                fprintf('%d,',i);
                [~, obj.hfile(i)] = ns_OpenFile([obj.recordingDir filesep obj.dataFileNames{i}]); %opens the mcd file for reading, which is used to avoid folder error in unix - bottle neck for this process
                [~, obj.fileInfo{i}] = ns_GetFileInfo(obj.hfile(i));
                
                obj.startDateLocal(i)=datenum(obj.fileInfo{i}.Time_Year,obj.fileInfo{i}.Time_Month,obj.fileInfo{i}.Time_Day,...
                    obj.fileInfo{i}.Time_Hour,obj.fileInfo{i}.Time_Min,obj.fileInfo{i}.Time_Sec+obj.fileInfo{i}.Time_MilliSec/1000);
                
                obj.recordingDurationLocal_ms(i)=obj.fileInfo{i}.TimeSpan*1000;
                obj.endDateLocal(i)=obj.startDateLocal(i)+datenum(0,0,0,0,0,obj.recordingDurationLocal_ms(i)/1000);
                
                [~, obj.entityInfo{i}] = ns_GetEntityInfo(obj.hfile(i), 1:obj.fileInfo{i}.EntityCount);
                
                obj.entityCount{i} = {obj.entityInfo{i}.ItemCount};
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
            
            obj.recordingDuration_ms=sum(obj.recordingDurationLocal_ms);
            obj.cumEnd=cumsum(obj.recordingDurationLocal_ms);
            obj.cumStart=[0 obj.cumEnd(1:end-1)];
            obj.startDate=obj.startDateLocal(1);
            obj.endDate=obj.endDateLocal(end);
            
            if obj.nRecordings>1
                obj.multifileMode=true;
            else
                obj.multifileMode=false;
            end
            
            % These properties are considered to be the same for all recorded files and extracted only from the first file
            obj.samplingFrequency=1/obj.fileInfo{1}.TimeStampResolution; %[Hz]
            obj.samplesPerMS=obj.samplingFrequency/1000;
            obj.entityLabels = {obj.entityInfo{1}.EntityLabel};
            obj.entityTypes = {obj.entityInfo{1}.EntityType};
            obj.entityName = cellfun(@(x) x(1:8),obj.entityLabels,'UniformOutput',0);
            
            obj.rawElectrodeEnteties = find(strcmp(obj.entityName,obj.defaultRawDataName));
            if isempty(obj.rawElectrodeEnteties)
                disp('Warming!!!! No raw electrode entities in recording');
                obj.rawElectrodeInfo = [];
            else
                obj.channelNames = strtrim(cellfun(@(x) x(25:27),obj.entityLabels(obj.rawElectrodeEnteties),'UniformOutput',0));
                %I changed this to get it to work!!!!
                
                % obj.channelNumbers = 1+str2double(cellfun(@(x) x(10:13),obj.entityLabels(obj.rawElectrodeEnteties),'UniformOutput',0));
                obj.channelNumbers=1:numel(obj.rawElectrodeEnteties);
                
                [ns,obj.rawElectrodeInfo] = ns_GetAnalogInfo(obj.hfile(1),obj.rawElectrodeEnteties);
                obj.rawChNumber2EntityNumber=zeros(1,max(obj.channelNumbers));
                obj.rawChNumber2EntityNumber(obj.channelNumbers)=obj.rawElectrodeEnteties;
                
                %translation between ch numbers and entities
                %!!! one is added to avoid probalems arrising with channel 0 - this is also reflected in channel maps
            end
            
            obj.ZeroADValue=-obj.rawElectrodeInfo(1).MinVal/obj.rawElectrodeInfo(1).Resolution;
            obj.dynamicRange=(obj.rawElectrodeInfo(1).MaxVal-obj.rawElectrodeInfo(1).MinVal)/obj.rawElectrodeInfo(1).Resolution;
            obj.MicrovoltsPerAD=obj.rawElectrodeInfo(1).Resolution*1e6;
            
            obj.trigggerEnteties = find(strcmp(obj.entityName,obj.defaultTriggerName));
            if isempty(obj.trigggerEnteties)
                disp('Warming!!!! No trigger entities in recording');
                obj.triggerInfo = [];
            else
                obj.triggerNames = strtrim(cellfun(@(x) x(25:27),obj.entityLabels(obj.trigggerEnteties),'UniformOutput',0));
                [ns,obj.triggerInfo] = ns_GetEventInfo(obj.hfile(1),obj.trigggerEnteties);
            end
            
            obj.digitalDataEnteties = find(strcmp(obj.entityName,obj.defaultDigitalDataName));
            if isempty(obj.digitalDataEnteties)
                disp('Warming!!!! No digital electrode entities in recording');
                obj.digitalDataEnteties=[];
            else
                [ns,obj.digitalDataInfo] = ns_GetAnalogInfo(obj.hfile(1),obj.digitalDataEnteties);
                pValid=find(~isempty({obj.digitalDataInfo.MaxVal}));
                obj.digitalDataChannels = strtrim(cellfun(@(x) x(25:27),obj.entityLabels(obj.digitalDataEnteties(pValid)),'UniformOutput',0));
                obj.digitalDataEnteties=obj.digitalDataEnteties(pValid);
                obj.digitalDataInfo=obj.digitalDataInfo(pValid);
            end
            
            obj.totalChannels=numel(obj.rawElectrodeEnteties);
            if max(obj.channelNumbers)~=numel(obj.channelNumbers) || min(obj.channelNumbers)~=1
                obj.channelNumbers=1:numel(obj.channelNumbers);
            end
            
            if obj.totalChannels<=32
                load('layout_300_6x6_FlexMEA.mat');
            elseif obj.totalChannels<=60
                load('layout_200_8x8.mat');
            elseif obj.totalChannels<=120
                load('layout_200_12x12.mat');
            elseif obj.totalChannels<=252
                load('layout_100_16x16.mat');
            end
            obj.chLayoutNumbers=En;
            obj.chLayoutNames=Ena;
            
            %remove from layout channels not in the channel list
            unrecordedChannels=setdiff(En(~isnan(En)),obj.channelNumbers);
            for i=1:numel(unrecordedChannels)
                p=find(En==unrecordedChannels(i));
                obj.chLayoutNumbers(p)=NaN;
                obj.chLayoutNames{p}=[];
            end
                        
            %The file should be closed, and not closing it may cause problems if many such files are opended.
            %However, if closed, reading consecutively from the file takes a lot of time since the file has to be reopened each time.
            %consequently the file is closed in the delete method overloaded on the dataRecording object
        end
    end
end