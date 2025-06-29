classdef matlabVariableRecording < dataRecording
%this class has not been fully written, the problem is that the different data types should first be define.
%My conclusion is that it may be better to add an option to keep some data in memory within each class itself.

    properties
        analogData
        electrodeData
        digitalData
    end

    properties (SetAccess=protected) %these properties are not considered as meta data and have to be loaded each time

    end

    properties (Constant)
        signalBits % the quantization of the sampling card
        numberOfCharFromEndToBaseName = "" %not relevant for this class but is a requirement of the superclass
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

        function obj=extractMetaData(obj)
        end
    end

    methods (Hidden)
        %class constructor
        function obj = matlabVariableRecording(dataVariable, dataType, channelNumbers, sourceDataObj)
            %(A, "analogData", 1, NP)
            %Usage: obj = matlabVariableRecording(recordingFile)
            
            props=?dataRecording;
            pNonConstant=cell2mat({props.PropertyList.Constant});
            for i=find(~pNonConstant)
                obj.(props.PropertyList(i).Name)=sourceDataObj.(props.PropertyList(i).Name);
            end

            obj.(dataType)=dataVariable;

        end

    end

end
