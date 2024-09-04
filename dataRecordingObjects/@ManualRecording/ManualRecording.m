classdef ManualRecording < dataRecording
    properties
       
    end
    
    properties (SetAccess=protected) %these properties are not considered as meta data and have to be loaded each time

    end
    
    properties (Constant)

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
        end
        

        function obj=closeOpenFiles(obj) %clear all open file handles
           
        end
        
        function obj=extractMetaData(obj)
            
        end
        
     
        
        function obj=getRecordingFiles(obj,recordingFile,fileExtension)
        end

    end
    
    methods (Hidden)
        %class constructor
        function obj = ManualRecording(recordingFile)

            scrsz = get(0,'ScreenSize');
            hFigure = figure('Position',[scrsz(3)*0.01 scrsz(4)*0.07 scrsz(3)*0.98 scrsz(4)*0.85], ...
                'Name','Manual data recording', 'NumberTitle','off', 'MenuBar','none', 'Toolbar','none', 'HandleVisibility','off');

            mainGrid=uix.Grid('Parent', hFigure, 'Padding', 2, 'Spacing', 5);
            getDataVariablePush=uicontrol('Parent', mainGrid, 'Callback',{@CallbackgetDataVariablePush}, 'Style','push', 'String','getData - desktop var');
            getDataVariableEdit=uicontrol('Parent', mainGrid, 'Style','edit', 'String','');

            samplingRateText=uicontrol('Parent', mainGrid, 'Style','text', 'String','');
            samplingRateEdit=uicontrol('Parent', mainGrid, 'Style','edit', 'String','');

            


            getDataVariablePush=uicontrol('Parent', mainGrid, 'Callback',{@CallbackgetDataVariablePush}, 'Style','push', 'String','getData - desktop var');
            
            getTriggerVariablePush=uicontrol('Parent', mainGrid, 'Callback',{@CallbackRestartPlotPush}, 'Style','push', 'String','Restart plot');
            getTriggerRecClassPush=uicontrol('Parent', mainGrid, 'Callback',{@CallbackExportPlotPush}, 'Style','push', 'String','Export plot');

            set(AVG.hGen.generalGrid, 'Widths',-1,'Heights', -1 );

            %Usage: obj = NeuraLynxRecording(recordingFile)
            obj.datatype='int16';

            function CallbackPaddingSamplesEdit(hObj,event)
                AVG.Params.filterPadding=str2num(get(hObj,'string'));
                AVG.Params.paddingSamples=AVG.Params.filterPadding/1000*AVG.Params.samplingFrequency(1);
            end

        end

    end

end
