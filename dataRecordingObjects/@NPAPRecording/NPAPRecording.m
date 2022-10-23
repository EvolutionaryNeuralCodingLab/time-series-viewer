classdef NPAPRecording < dataRecording
    properties
        %{
        recordingName %(String) The name of the recording
        recordingDir % (String) Full directory containing the recorded session
        dataFileNames % (Cell 1 x N)  array of N recording data file names
        startDate %(1x1) Start date (time) of Recording (matlab date number format)
        endDate %(1x1) End date (time) of Recording (matlab long format)
        samplingFrequency %(1xN) Sampling rate [Hz]
        recordingDuration_ms %(1x1) the total duration of the recording in [ms]
        channelNames % (Cell 1xN) a cell array with the N names of the channels
        channelNumbers % (1xN) an array with integer channel numbers (>=1 integer)
        
        triggerNames %the names of trigger channels (not critical)
        analogChannelNumbers % (1xN) the numbers of channels containing non-electrode analog inputs
        analogChannelNames % (1xN) the numbers of channels containing non-electrode analog inputs
        
        chLayoutNumbers %(MxN) The layout of the channel numbers in physical space arranged in an M by N grid
        chLayoutNames %(Cell MxN)The layout of the channel names in physical space arranged in an M by N grid
        electrodePitch % distance between electrodes (not critical)
        chLayoutPositions % (1xN or 2xN or 3xN) array of electrode position in [x or x,y or x,y,z]
        layoutName %the name of the channel layout (electrode type)
        
        convertData2Double = 1; % if data should be converted to double from the original quantization
        ZeroADValue % the digital zero value
        MicrovoltsPerAD % the digital to analog conversion value
        datatype        % class of data in the recording
        %}
    end
    
    properties (SetAccess=protected) %these properties are not considered as meta data and have to be loaded each time
        
    end
    
    properties (Constant)
        defaultLocalDir='C:\Users\Tulip\Documents\Academic\Post-Doc\Experiments'; %Default directory from which search starts
    end

    methods
        
        function [V_uV,t_ms]=getData(obj,channels,startTime_ms,window_ms)
            %Extract Neuralynx recording data from file to memory
            %Usage: [V_uV,t_ms]=obj.getData(channels,startTime_ms,window_ms);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        t_ms - A time vector relative to recording start (t=0 at start)
            
            
        end
        function [V_uV,t_ms]=getAnalogData(obj,channels,startTime_ms,window_ms)
            %Extract Neuralynx recording data from file to memory
            %Usage: [V_uV,t_ms]=obj.getData(channels,startTime_ms,window_ms);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        t_ms - A time vector relative to recording start (t=0 at start)
        
        end
        
        function [T_ms,chNumber]=getTrigger(obj,startTime_ms,window_ms,name)
            %Extract triggers from file Neuralynx recording
            %Usage : [T_ms]=obj.getTrigger(startTime_ms,window_ms)
            %Input : name - which bit to extract for time stamps (out of 8,default = first bit, 1)
            %        startTime_ms - start time [ms].
            %        window_ms - the window duration [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %
            %Output: T_ms - trigger times [ms] - different triggers are arranged in a cell array
            
            %read Data from *.nev files (NeuroLynx) - files containing the data
            %fieldSelectionTrigger: (1): Timestamps, (2): Event IDs, (3): TTLs, (4): Extras, (5): Event Strings
            %Bits are coded as decimal integers
            
        end
        
        function obj=closeOpenFiles(obj) %clear all open file handles
            if ~iscell(obj.recordingDir)
                nRecordings=1;
            else
                nRecordings=numel(obj.recordingDir);
            end
            if ~isempty(obj.fid)
                for j=1:nRecordings
                    fclose(obj.fid(j));
                end
            end
            if ~isempty(obj.fidA)
                for j=1:nRecordings
                    fclose(obj.fidA(j));
                end
            end
            if ~isempty(obj.fidEvnt)
                fclose(obj.fidEvnt);
            end
        end
        
        function obj=extractMetaData(obj)
            fprintf('Extracting meta data from: %s...\n',obj.recordingDir);
            disp('saving meta data');
            obj.saveMetaData;
        end
        
    end
    
    methods (Hidden)
        %class constructor
        function obj = NPAPRecording(recordingFolder)
            %Usage: obj = NeuraLynxRecording(recordingFile)
            obj.folderMode=1;
            %get data files
            if nargin==0
                recordingFolder=[];
            elseif nargin>1
                disp('Object was not constructed since too many parameters were given at construction');
                return;
            end
            %obj=obj.getRecordingFiles(recordingFolder);
                        
            if exist(obj.metaDataFile,'file') && ~obj.overwriteMetaData
                obj=loadMetaData(obj);
            else
                obj=extractMetaData(obj);%file identifiers are aquired inside function
            end
            
            %load layout
            obj=obj.loadChLayout;
            if isempty(obj.chLayoutNumbers)
                disp('No layout file in directory, to add layout please add a *.chMap file to the data directory with the layout name');
            end
            
        end
        
    end
    
end

%}