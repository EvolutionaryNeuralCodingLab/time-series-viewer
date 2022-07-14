classdef HD3BrainRecording < dataRecording
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
    end
    properties (Constant, Hidden)
        fileExtension='brw';
        defaultLocalDir='C:\Users\Tulip\Documents\Academic\Post-Doc\Experiments'; %Default directory from which search starts
    end
    properties (Hidden)
        Ch
        nCh
        NFrames
        NStreams
        MicrovoltsPerAD
        recordingDurationLocal_ms
        brwRdr
        ExistRawStream
        ExistSpikeStream
    end
    methods
        function [V_uV,T_ms]=getData(obj,channels,startTime_ms,window_ms,name)
            %Extract MCRack recording data from file to memory
            %Usage: [V_uV,T_ms]=obj.getData(channels,startTime_ms,window_ms,name);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the stream (if not entered, default name is used)
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        T_ms - A time vector relative to recording start (t=0 at start)
            if nargin<4
                error('method getData was not used correctly: wrong number of inputs');
            end
            if ~obj.ExistRawStream
                error('No raw data stream exists in this recording');
            end
            
            conversionFactor=1/1000*obj.samplingFrequency;
            startSample=round(startTime_ms*conversionFactor);
            startTime_ms=startSample/conversionFactor;
            
            windowSamples=round(window_ms*conversionFactor);
            window_ms=windowSamples/conversionFactor;
            
            endTime_ms=startTime_ms+window_ms; %no need to conversion factor
            recordingDuration_ms=round(obj.recordingDuration_ms*conversionFactor)/conversionFactor;
           
            nChannels=numel(channels);
            nTrials=length(startTime_ms);
            V_uV=zeros(nChannels,nTrials,windowSamples);
            
            %cumStart=[-Inf obj.cumStart obj.cumEnd(end)];
            %cumEnd=[0 obj.cumEnd Inf];
            
            
            if obj.multifileMode %this mode currently does not support extraction from edges of the recording
               
            else
                recFile=1;
                for i=1:nTrials
                    if startTime_ms(i)>=0 && (startTime_ms(i)+window_ms)<=recordingDuration_ms
                        if nChannels==obj.nCh % for large numbers of channels-> reads all channels first and then deletes unwanted ones (50 was chosen since the speed for this was optimal
                            V_uV(:,i,:) = double(obj.brwRdr{recFile}.GetRawData(obj.Ch{recFile}, startSample(i), windowSamples));% load raw data
                        else
                            chsSubset = NET.createArray('BW.ChCoord', nChannels);
                            for j=1:nChannels
                                chsSubset(j)=obj.Ch{recFile}(channels(j));
                            end
                            V_uV(:,i,:) = double(obj.brwRdr{recFile}.GetRawData(chsSubset, startSample(i), windowSamples));% load raw data
                        end
                    else
                        startSample=min(0,round(startTime_ms(i)*conversionFactor));
                        endSample=min(windowSamples,round((recordingDuration_ms-startTime_ms(i))*conversionFactor)); %in window (not in recording)
                        if nChannels>50 % for large numbers of channels-> reads all channels first and then deletes unwanted ones (50 was chosen since the speed for this was optimal
                            tmpData = double(obj.brwRdr{recFile}.GetRawData(obj.Ch{recFile}, startSample, endSample));% load raw data
                            tmpData = tmpData(channels,:);
                        else
                            for j=1:nChannels
                                tmpData = double(obj.brwRdr{recFile}.GetRawData(obj.Ch{recFile}(channels(j)), startSample, endSample));% load raw data
                            end
                        end
                        V_uV(:,i,1-startSample:endSample)=tmpData;
                        disp('Recording at edge');
                    end
                end
            end
            
            if nargout==2
                T_ms=(1:windowSamples)*(1e3/obj.samplingFrequency);
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
            D=[];
            T_ms={};
        end
        function [T_ms]=getTrigger(obj,startTime_ms,window_ms,name)
            %Extract triggers from file Neuralynx recording
            %Usage : [T_ms]=obj.getTrigger(,startTime_ms,window_ms,name)
            %Input : startTime_ms - start time [ms].
            %        window_ms - the window duration [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        name
            %Output: T_ms - trigger times [ms] - different triggers are arranged in a cell array
            T_ms={};
        end
        
        function obj=closeOpenFiles(obj)
            for i=1:obj.nRecordings
                obj.brwRdr{i}.Close;
            end
        end
        
    end %Methods

    methods (Hidden)
        %class constructor
        function obj = HD3BrainRecording(recordingFile)
            %Usage: obj = MCRackRecording(recordingFile),
            %structure of resulting file: [nChannels x nTrials x nSamples]
            
            %get data files
            if nargin==0
                recordingFile=[];
            elseif nargin>1
                disp('Object was not constructed since too many parameters were given at construction');
                return;
            end
            obj=obj.getRecordingFiles(recordingFile,obj.fileExtension);
            
            libPath = which('BrwExtReader.dll'); %find library of brw file reader
            asm = NET.addAssembly(libPath); % load the private assembly by the full path (comprising file name and extension)
            try
                %open files for reading and getting file data
                for i=1:obj.nRecordings
                    obj.brwRdr{i}=BW.BrwRdr; % create an instance of the BrwFile class
                    status = obj.brwRdr{i}.Open(obj.recordingName{i});
                    
                    obj.recordingDurationLocal_ms(i) = obj.brwRdr{i}.RecDuration;
                    obj.NFrames(i)=obj.brwRdr{i}.RecNFrames;
                    
                    obj.Ch{i}=obj.brwRdr{i}.GetRecChsUnion();%  gets the union of all the recorded channels regardless the recorded stream
                end
                %methodsview(obj.brwRdr{i})
                obj.recordingDuration_ms=sum(obj.recordingDurationLocal_ms);
                
                obj.samplingFrequency=double(obj.brwRdr{1}.SamplingRate); %assuming that sampling rate did not change across sessions
                obj.ExistRawStream=obj.brwRdr{1}.ExistRawStream;
                obj.ExistSpikeStream=obj.brwRdr{1}.ExistSpikeStream;
                
                if obj.nRecordings>1
                    obj.multifileMode=true;
                else
                    obj.multifileMode=false;
                end

                obj.nCh = double(obj.Ch{1}.Length);
                
                obj.channelNumbers=1:obj.nCh;
                obj.channelNames=cell(1,obj.nCh);
                rows=zeros(1,obj.nCh);
                cols=zeros(1,obj.nCh);
                
                for i = 1:obj.nCh
                    obj.channelNames{i} = char(obj.Ch{1}(i).ToString());
                    rows(i)=double(obj.Ch{1}(i).Row);
                    cols(i)=double(obj.Ch{1}(i).Col);
                end
                NRows=obj.brwRdr{1}.MeaDeviceLayoutNRows;
                NColumns=obj.brwRdr{1}.MeaDeviceLayoutNColumns;
                
                obj.chLayoutNumbers=zeros(NRows,NColumns);
                obj.chLayoutNames=cell(NRows,NColumns);
                
                obj.chLayoutNumbers(sub2ind([NRows NColumns],rows,cols))=obj.channelNumbers;
                obj.chLayoutNames(sub2ind([NRows NColumns],rows,cols))=obj.channelNames;

                %obj.ZeroADValue
            catch Err
                for i=1:obj.nRecordings
                    obj.brwRdr{i}.Close;
                end
                disp({'An error occured while creating Brain Wave Reader object';Err.message});
            end
        end
    end
end