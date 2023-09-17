classdef AVPlotData < handle
    properties
        layoutName={};          %layout of electrode if existing 
        T=[];                   %time vector
        startTime=[];           %recording start time
        nCh=[];                 %number of channels
        nTrials=[];             %number of trials
        nSamples=[];            %number of samples
        channelNames=[];        %channel names
        channelNumbersOld=[];   %old channel numbers - replica of channel number that allways keeps the previous value
        hPlotAxis='empty';           %handle for the axis for plotting
        hControlPanel='empty';       %handle for the panel on which controls will be created
        hPlot=[];               %handle for the plot data
        hPlotControls=[];       %handle for plot controls
        hPlotBackground=[];     %handle for plot background (not deleted when plotting old background on new plot
        plotParams=[];          %structure with plot specific parameters
        currentPlotName=[];     %the name of the current from the plot list (gathered with method getPlotNames)
        plotMethods=[];         %plot methods extracted from plot library
        plotCreateMethods=[];   %plot create methods
        plotNames=[];           %plot names extracted from plot library
        currentPlotIdx=[];
        refreshPlot=0;
        plotAnalogChannels=0;   %whether to plot the analog channels
        analogChannelNames=[];  %analog channel names
        analogScaling=1;
        yl=[];                  %the data limits used in plots
    end
    properties (Hidden)
        recordingObjReference=[]; %a handle to a recording object to modify recording properties directly
    end
    properties (SetObservable)
        M=[];                   %activity data [nChannels x nTrials x nSamples]
        channelNumbers=[];      %channel numbers
        A=[];                   %analog data [nChannels x nTrials x nSamples]
    end
    methods
        %class constractor
        function obj=AVPlotData(startMethod)
            addlistener(obj,'M','PostSet',@obj.changedActivityData); %add a listener to M, after its changed its size is updated in the changedDataEvent method
            addlistener(obj,'channelNumbers','PostSet',@obj.changedChannelNumber); %add a listener to M, after its changed its size is updated in the changedDataEvent method
            %get the plot names
            allMethods=methods(class(obj));
            obj.plotMethods=allMethods(strncmp('plot',allMethods,4));
            obj.plotNames=cellfun(@(x) x(5:end),obj.plotMethods,'UniformOutput',0); %remove the "plot" from the string
            obj.plotCreateMethods=cellfun(@(x) ['create' x(5:end)], obj.plotMethods,'UniformOutput',0);
            if nargin>=1
                obj.(startMethod);
            end
        end
        function GUI(obj)
            % have to write the code - most lines below should be revised
            
            scrsz = get(0,'ScreenSize');
            hFigure=figure('Position',[scrsz(3)*0.01 scrsz(4)*0.07 scrsz(3)*0.98 scrsz(4)*0.85], 'Name','AVPlotData GUI', 'NumberTitle','off');%,'CloseRequestFcn',@closeMainGUIFigure);
            
            % set file menus
            hFile = uimenu(hFigure, 'Label', 'File' );
            AVG.hMenubar.exportData=uimenu( AVG.hMenubar.hFile, 'Label', 'Export data', 'Callback', @CallbackExportData);
            %AVG.hMenubar.fileLoad=uimenu( AVG.hMenubar.hFile, 'Label', 'Load', 'Callback', @CallbackMenuLoad);
            %AVG.hMenubar.fileExit=uimenu( AVG.hMenubar.hFile, 'Label', 'Exit', 'Callback', @CallbackMenuExit);
            
            % set plot meuns
            AVG.hMenubar.hPlotType = uimenu(AVG.hMainFigure.hFigure, 'Label', 'Plot' );
            for i=1:length(AVG.Params.PlotList)
                AVG.hMenubar.hPlot.(AVG.Params.PlotList{i})=uimenu(AVG.hMenubar.hPlotType,...
                    'Label', AVG.Params.PlotList{i}, 'Checked','off', 'Callback', @CallbackPlotTypeSelection);
            end
            set(AVG.hMenubar.hPlot.(AVG.Params.currentPlotName), 'Checked','on');
            
            initializePlot(obj,plotName,hPlotAxis,hControlPanel);
            
        end
        %all class methods are written as separate functions
        function [plotMethods,plotNames]=getPlotNames(obj)
            plotMethods=obj.plotMethods;
            plotNames=obj.plotNames;
        end
        function changedActivityData(obj,event,metaProp) %refresh properties every time new data are loaded
            oldnTrials=obj.nTrials;
            [obj.nCh obj.nTrials obj.nSamples]=size(obj.M);
            if oldnTrials~=obj.nTrials %a change in the number of trials requires replotting (this does not hold for change in the time stamp value) 
                obj.refreshPlot=1;
            end
        end
        function changedChannelNumber(obj,event,metaProp) %define changes to plot every time the channel configuration is changed
            if ~isequal(obj.channelNumbers,obj.channelNumbersOld)
                obj.refreshPlot=1;
            end
            obj.channelNumbersOld=obj.channelNumbers;
        end
        function deletePlotData(obj)
            delete(obj.hPlot);
            obj.hPlot=[];
            set(obj.hPlotAxis,'YDir','normal');
        end
        function deletePlotControls(obj)
            if isfield(obj.hPlotControls,'plotPropGrid')
                delete(obj.hPlotControls.plotPropGrid);
                obj.hPlotControls=[]; %notice that after deleting plotPropGrid all controls are also deleted (children)
            end
        end
        function deletePlotBackground(obj)
            delete(findobj('Childred',obj.hPlotAxis));
            %{
            if ~isempty(obj.hPlotBackground)
                delete(cell2mat(struct2cell(obj.hPlotBackground)));
                obj.hPlotBackground=[];
            end
            %}
        end
        function generateFalsePlotData(obj)
            obj.hPlot=plot(obj.hPlotAxis,1,1);
        end
        function initializePlot(obj,plotName,hPlotAxis,hControlPanel)
            %set handles to main plot and plot controls
            if nargin==2 %if evoked independent of activity viewer
                fControlPanel = figure('Position',[100 100 400 700],'name','AVDataPlot - controls');
                hControlPanel = uix.Panel('Parent',fControlPanel, 'Title','Plot specific properties');
                fPlotAxis=figure('name','AVDataPlot - plot','position',[530 100 1000 700]);
                hPlotAxis = axes('Parent',fPlotAxis);
            elseif nargin==3
                fPlotAxis=figure('name','AVDataPlot - plot','position',[530 100 1000 700]);
                hPlotAxis = axes('Parent',fPlotAxis);
            end
            obj.plotParams=[];
            obj.hPlotAxis=hPlotAxis;
            set(obj.hPlotAxis,'XTickMode','auto','YTickMode','auto','yscale','lin');
            obj.hControlPanel=hControlPanel;
            obj.currentPlotIdx=strcmp(obj.plotNames,plotName);
            obj.(obj.plotCreateMethods{obj.currentPlotIdx});
        end
        function replot(obj)
            obj.deletePlotData;
            obj.(obj.plotMethods{obj.currentPlotIdx});
        end
    end
end %EOF