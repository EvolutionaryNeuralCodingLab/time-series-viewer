classdef AIMatlab < dataRecording
    properties
        pMultiSeries=[];
        t %ms
        ic
        I
        trig %[ms]
        estimatedResolution %ms
        maxSamplingFrequency
    end
    properties (Constant, Hidden)
        fileExtension='mat';
    end
    properties (Constant)
        defaultLocalDir='C:\Users\Tulip\Documents\Academic\Post-Doc\Experiments'; %Default directory from which search starts
        signalBits                   	% defined in dataRecording
        numberOfCharFromEndToBaseName	% defined in dataRecording
    end
    properties (Hidden)
        
    end
    properties (SetObservable)
        timeBin = 1; %ms - defines the sampling frequency
    end
    methods
        function reloadMetaData(obj)
            obj.recordingDuration_ms=max(obj.t);
            nNeurons=size(obj.ic,2);
            for i=1:nNeurons
                obj.channelNames{i}=[num2str(obj.ic(1,i)) '-' num2str(obj.ic(2,i))];
            end
            obj.channelNumbers=1:nNeurons;
        end
        function obj=updateSamplingFrequency(obj,event,metaProp)
            obj.samplingFrequency=1/obj.timeBin*1000;
        end
        function obj=updateMaximalSamplingFrequency(obj)
            if obj.dataFileNames{1}~=0
                load([obj.recordingDir filesep obj.dataFileNames{1}],'-mat', 't');
            else
                t=obj.t;
            end
            if ~isempty(obj.pMultiSeries)
                t=t{obj.pMultiSeries};
            end
            %estimate resolution from recorind
            tmp=diff(sort(t));
            obj.estimatedResolution=min(abs(tmp(tmp>0)));
            obj.maxSamplingFrequency=1000/obj.estimatedResolution; %look for the minimal time difference between samples that is larger than 0
            if obj.maxSamplingFrequency>10000 %for cases where there are rounding issues
                obj.maxSamplingFrequency=round(obj.maxSamplingFrequency); 
                obj.estimatedResolution=1000/obj.maxSamplingFrequency;
            end
            obj.t=t;
        end
        function [V_uV,t_ms]=getData(obj,channels,startTime_ms,window_ms)
            %Extract filtered matlab data from file to memory
            %Usage: [V_uV,t_ms]=obj.getData(channels,startTime_ms,window_ms);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials  
            %        t_ms - A time vector relative to recording start (t=0 at start)
            [V_uV]=permute(BuildBurstMatrixA(obj.ic(:,channels),round(obj.t/obj.timeBin),obj.I,round(startTime_ms/obj.timeBin),round(window_ms/obj.timeBin)),[2 1 3]);
            nSamples=size(V_uV,3);
            if nargout==2
                t_ms=(1:nSamples)*obj.timeBin;
            end
        end
        function [T_ms]=getTrigger(obj,startTime_ms,window_ms,name)
            T_ms={};
            load([obj.recordingDir obj.dataFileNames],'trigger');
            if exist('trigger','var')
                if ~iscell(trigger)
                    T_ms{1}=trigger;
                else
                    T_ms=trigger;
                end
            end
        end
    end
    
    methods (Hidden)
        %class constructor
        function obj = AIMatlab(recordingFile)
            
            %get data files
            if nargin==0
                recordingFile=[];
            elseif nargin>1
                disp('Object was not constructed since too many parameters were given at construction');
                return;
            end
            obj=obj.getRecordingFiles(recordingFile,obj.fileExtension);
            if obj.dataFileNames{1}==0
                tvar = uigetvariables({'Select time variable','Select index channel variable'});
                obj.t=tvar{1};
                obj.ic=tvar{2};
                ic=obj.ic;
            else
                load([obj.recordingDir filesep obj.dataFileNames{1}],'ic','I'); %loads t

            end
            
            if iscell(ic) %in case t/I/ic are structured as cell array containning more than one data set
                nSeries=numel(ic);
                f=figure('position',[100 100 100 500],'name','Select series');
                VB = uiextras.VBox('Parent',f, 'Spacing',4, 'Padding',4);
                uicontrol('parent',VB,'style','text','string',{'File contains more than one data sets','Please select dataset number'},'horizontalAlignment','left');
                h=uicontrol('parent',VB,'style','listbox','string',1:nSeries);
                set(VB,'Sizes',[-1 -10]);
                %waitfor(h,'value'); % does not continue if the default option was pressed
                waitforbuttonpress;
                pause(0.1); %wait for handle h to update, otherwise value sometimes is not updated withn until next line
                
                obj.pMultiSeries=get(h,'value');
                delete(f);
                ic=ic{obj.pMultiSeries}; % t is loaded later
                if exist('I','var')
                    I=I{obj.pMultiSeries};
                end
            end
            
            addlistener(obj,'timeBin','PostSet',@obj.updateSamplingFrequency);
            
            obj=obj.updateMaximalSamplingFrequency;
            obj.samplingFrequency=1/obj.timeBin*1000;
            obj.ic=ic;
            if exist('I','var')
                obj.I=I;
            else
                obj.I=ones(1,numel(obj.t));
            end
            reloadMetaData(obj);
        end
    end
end