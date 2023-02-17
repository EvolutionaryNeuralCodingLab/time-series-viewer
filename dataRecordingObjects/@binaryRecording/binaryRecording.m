classdef binaryRecording < dataRecording
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
        bytesPerSample
        nTotSamples
        totalChannels
        totalAnalogChannels
    end
    
    properties (Constant)
        defaultLocalDir='C:\Users\Tulip\Documents\Academic\Post-Doc\Experiments'; %Default directory from which search starts
        signalBits = 16; %the quantization of the sampling card
        fileExtension='bin';
        numberOfCharFromEndToBaseName %the number of characters from the end of the file name to get to the base file name
    end
    
    properties (Hidden)
        fid
        fidAnalog
    end
    
    methods
        
        function [V_uV,T_ms]=getData(obj,channels,startTime_ms,window_ms,name)
            %Extract MCRack recording raw data from file
            %Usage: [V_uV,T_ms]=obj.getData(channels,startTime_ms,window_ms,name);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the stream (if not entered, default name is used)
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        T_ms - A time vector relative to recording start (t=0 at start)
            if nargin<4
                error('wrong number of inputs,use: [V_uV,T_ms]=getData(obj,channels,startTime_ms,window_ms)');
            end
            if isempty(channels) %if no channels are entered, get all channels
                channels=obj.channelNumbers;
            end
            conversionFactor=1/1000*obj.samplingFrequency;
            startTime_ms=round(startTime_ms*conversionFactor)/conversionFactor;
            window_ms=round(window_ms*conversionFactor)/conversionFactor;
            endTime_ms=startTime_ms+window_ms; %no need to conversion factor
            windowSamples=round(window_ms*conversionFactor);
            startSample=round(startTime_ms*conversionFactor);
            nTrials=length(startTime_ms);
            V_uV=ones(numel(channels),nTrials,windowSamples,obj.datatype);
                        
            for i=1:nTrials
                if startTime_ms(i)>=0 && (startTime_ms(i)+window_ms)<=obj.recordingDuration_ms
                    fseek(obj.fid,startSample(i)*obj.bytesPerSample*obj.totalChannels,'bof');
                    data = fread(obj.fid,obj.totalChannels*windowSamples,['*' obj.datatype]);
                    data = reshape(data, obj.totalChannels,windowSamples);
                    V_uV(:,i,:)=data(channels,:);
                else
                    startSampleTmp=min(0,startSample(i));
                    endSampleTmp=min(windowSamples,obj.nTotSamples-startSample(i)); %end sample in window (not in recroding)
                    
                    startSampleRec=max(0,startSample(i));
                    endSampleRec=min(startSample(i)+windowSamples,obj.nTotSamples);
                    
                    fseek(obj.fid,startSampleRec*obj.bytesPerSample*obj.totalChannels,'bof');
                    data = fread(obj.fid,(endSampleRec-startSampleRec)*obj.totalChannels,['*' obj.datatype]);
                    data = reshape(data, obj.totalChannels,endSampleRec-startSampleRec);
                    
                    V_uV(:,i,1-startSampleTmp:endSampleTmp)=data(channels,:);
                    disp('Recording at edge');
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
            if nargin<4
                error('wrong number of inputs,use: [V_uV,T_ms]=getAnalogData(obj,channels,startTime_ms,window_ms)');
            end
            if isempty(channels) %if no channels are entered, get all channels
                channels=obj.analogChannelNumbers;
            end
            conversionFactor=1/1000*obj.samplingFrequencyAnalog;
            startTime_ms=round(startTime_ms*conversionFactor)/conversionFactor;
            window_ms=round(window_ms*conversionFactor)/conversionFactor;
            endTime_ms=startTime_ms+window_ms; %no need to conversion factor
            windowSamples=round(window_ms*conversionFactor);
            startSample=round(startTime_ms*conversionFactor);
            nTrials=length(startTime_ms);
            V_uV=ones(numel(channels),nTrials,windowSamples,'int16');
                        
            for i=1:nTrials
                if startTime_ms(i)>=0 && (startTime_ms(i)+window_ms)<=obj.recordingDuration_ms
                    fseek(obj.fidAnalog,startSample(i)*obj.bytesPerSample*obj.totalAnalogChannels,'bof');
                    data = fread(obj.fidAnalog,obj.totalAnalogChannels*windowSamples,['*' obj.datatype]);
                    data = reshape(data,obj.totalAnalogChannels,windowSamples);
                    V_uV(:,i,:)=data(channels,:);
                else
                    startSample=min(0,startSample(i));
                    endSample=min(windowSamples,obj.nTotSamples-startSample(i)); %end sample in window (not in recroding)
                    
                    startSampleRec=max(0,startSample(i));
                    endSampleRec=min(startSample(i)+windowSamples,obj.nTotSamples);
                    
                    fseek(obj.fidAnalog,startSampleRec*obj.bytesPerSample*obj.totalAnalogChannels,'bof');
                    data = fread(obj.fid,(endSampleRec-startSampleRec)*obj.totalAnalogChannels,['*' obj.datatype]);
                    data = reshape(data(channels,:), obj.totalAnalogChannels,endSample+startSample);
                    
                    V_uV(:,i,1-startSample:endSample)=data;
                    disp('Recording at edge');
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
            disp('Method not implemented yet!!!');
        end
        
        function [T_ms]=getTrigger(obj,startTime_ms,window_ms,name)
            %Extract triggers from recording
            %Usage : [T_ms]=obj.getTrigger(,startTime_ms,window_ms,name)
            %Input : startTime_ms - start time [ms].
            %        window_ms - the window duration [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        name
            %Output: T_ms - trigger times [ms] - different triggers are arranged in a cell array
            triggerFile=[obj.recordingDir filesep obj.recordingName '_Triggers.' obj.fileExtension];
            if exist(triggerFile,'file')
                fid=fopen(triggerFile,'r');
                nTriggers=double(fread(fid,numel(obj.triggerNames),'*uint32')-1);
                for i=1:numel(nTriggers)
                   	T_ms{i}=double(fread(fid,nTriggers(i),'*uint32')-1)'/(obj.samplingFrequency/1000);
                end
            end
            fclose(fid);
        end
        
        function obj=closeOpenFiles(obj)
            fclose(obj.fid);
            if obj.fidAnalog>0
                fclose(obj.fidAnalog);
            end
        end
        
        function [obj]=extractMetaData(obj)
            %get data from meta file (jrclust fields)
            %The following fields are required: 'sRateHz','scale','nChans'
            metaFileName=obj.recordingName;
            fullMetaFileName=[obj.recordingDir filesep metaFileName '_meta.txt'];
            if exist(fullMetaFileName,'file')
                fidMeta=fopen(fullMetaFileName,'r');
                metaData = textscan(fidMeta,'%s = %s','TextType','string','Delimiter',' = ');  %'Whitespace','',
                fclose(fidMeta);
                
                for i=1:size(metaData{1},1)
                    if any(isletter(metaData{2}{i}))
                        T.(metaData{1}{i})=metaData{2}{i};
                    else
                        T.(metaData{1}{i})=str2num(metaData{2}{i});
                    end
                end
                obj.samplingFrequency=T.sRateHz;
                obj.samplingFrequencyAnalog=T.sRateAnalogHz;
                obj.MicrovoltsPerAD=T.scale;
                obj.MicrovoltsPerADAnalog=T.scaleAnalog;
                obj.ZeroADValue=T.zeroADValue;
                obj.ZeroADValueAnalog=T.zeroADValueAnalog;
                obj.totalChannels=T.nChans;
                obj.datatype=T.vcDataType;
                obj.totalAnalogChannels=T.nAnalogChans;
                obj.triggerNames=mat2cell(1:T.nTriggerChans,1,ones(1,T.nTriggerChans));
                if isfield(T,'channelNumbersAnalog') %for legacy version when channelNumbersAnalog was not saved
                    obj.analogChannelNumbers=1:T.nAnalogChans;
                end
                obj.channelNumbers=T.channelNumbers;
                obj.channelNames=cellfun(@(x) num2str(x),mat2cell(obj.channelNumbers,1,ones(1,numel(obj.channelNumbers))),'UniformOutput',0); %
                obj.analogChannelNames=cellfun(@(x) num2str(x),mat2cell(obj.analogChannelNumbers,1,ones(1,numel(obj.analogChannelNumbers))),'UniformOutput',0);

                if isfield(T,'vcProbe')
                    %if numel(T.vcProbe)<7
                    %    error('Probe name should have the prefix "Layout_". Please correct and run again.');
                    %end
                    if strcmp(T.vcProbe(1:min(7,numel(T.vcProbe))),'Layout_')
                        obj=obj.loadChLayout(T.vcProbe(8:end));%remove the layout ending
                    else
                        obj=obj.loadChLayout(T.vcProbe);%remove the layout ending
                    end
                end
                %check the number of samples in the binary file
                fid=fopen([obj.recordingDir filesep obj.dataFileNames{1}],'r');
                fseek(fid, 0, 'eof');
                position = ftell(fid);
                obj.nTotSamples=floor(position/2/obj.totalChannels);
                fclose(fid);
                
                obj.recordingDuration_ms=obj.nTotSamples/obj.samplingFrequency(1)*1000;
                disp('saving meta data');
                obj.saveMetaData;
            else
                error('could not read meta data file (.meta) for recording, data object not created!!!');
            end
        end
    end
    
    methods (Hidden)
        %class constructor
        function obj = binaryRecording(recordingFile)
            %Usage: obj = MCRackRecording(recordingFile),
            %structure of resulting file: [nChannels x nTrials x nSamples]

            obj.datatype='int16';
            obj.bytesPerSample=2;
            obj.multifileMode=false;
            obj.ZeroADValue=0;
            obj.ZeroADValueAnalog=0;
            
            %get data files
            if nargin==0
                recordingFile=[];
            elseif nargin>1
                disp('Object was not constructed since too many parameters were given at construction');
                return;
            end
            if isstr(recordingFile)
                tmp=recordingFile;
                clear recordingFile;
                recordingFile{1}=tmp;
            end
            if ~any(strcmp(recordingFile{1}(end-2:end),{'bin','dat'}))
                warning('Recording file is not given, this may create errors, please provide a *.bin filename');
            end
            obj=obj.getRecordingFiles(recordingFile,obj.fileExtension);
            
            if isfile([obj.metaDataFile,'.mat']) && ~obj.overwriteMetaData
                obj=loadMetaData(obj);
            else
                obj=extractMetaData(obj);
            end
            
            obj=obj.loadChLayout;
            
            obj.fid=fopen([obj.recordingDir filesep obj.recordingName '.' obj.fileExtension],'r');
            if numel(obj.analogChannelNumbers)>0
                obj.fidAnalog=fopen([obj.recordingDir filesep obj.recordingName '_Analog.' obj.fileExtension],'r');
            end
        end
    end
end