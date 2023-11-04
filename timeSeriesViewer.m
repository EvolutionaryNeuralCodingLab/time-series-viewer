function timeSeriesViewer(varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%% Definition of default values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AVG.recordingObj=[]; %initiate recording object
AVG.hVideoSyncFigure.hFigure=[];
AVG.Params.timeWindowExternalHandles=[]; %to sync gui to another one
AVG.Params.videoSyncFile=[]; %to sync with video file
AVG.recordingObj=[]; %recording object
AVG.Params.videoReader=[]; %video reader for syncing video to ephys

%directories
AVG.Params.NSKToolBoxMainDir=fileparts(which('identifierOfMainDir4NSKToolBox'));
AVG.Params.timeSeriesViewerMainDir=fileparts(which('timeSeriesViewer.m'));

%Check that GUI layout is installed and if not install it
installGUILayoutToolBox=false;
toolboxes = matlab.addons.toolbox.installedToolboxes;
if isempty(toolboxes)
    installGUILayoutToolBox=1;
else
    toolboxesNames={toolboxes.Name};
    if ~any(strcmp(toolboxesNames,'GUI Layout Toolbox'))
        installGUILayoutToolBox=1;
    end
end
if installGUILayoutToolBox
    disp('GUI Layout toolbox is not installed, trying to install uix version 2.3.4...');
    d=which('GUI Layout Toolbox 2.3.4.mltbx');
    installedToolbox = matlab.addons.toolbox.installToolbox(d,true);
end

if isempty(ver('layout'))
    if ~downloadGUILayout
        return;
    end
end

dirSep=filesep; %choose file/dir separator according to platform

AVG.Params.dataClassDirectory=[AVG.Params.timeSeriesViewerMainDir dirSep 'dataRecordingObjects' dirSep];
AVG.Params.dataAnalysisDirectory=[AVG.Params.timeSeriesViewerMainDir dirSep 'dataAnalysis' dirSep];
AVG.Params.colormapFileDirectory=[AVG.Params.timeSeriesViewerMainDir dirSep 'colormaps' dirSep];

AVG.Params.defaultColorMap='colormapHot128';
colormapFiles=dir([AVG.Params.colormapFileDirectory 'colormap*.mat']);
AVG.Params.colormapList=cellfun(@(x) x(1:end-4),{colormapFiles.name},'UniformOutput',0);
AVG.Params.currentColormap=find(strcmp(AVG.Params.colormapList,AVG.Params.defaultColorMap));
AVG.Params.colormap=load(AVG.Params.colormapList{AVG.Params.currentColormap},'cmap');
AVG.Params.colormap=AVG.Params.colormap.cmap;

AVG.Params.defaultActiveChannels=[];
%AVG.Params.defaultPlotName='ActivityTracePhysicalSpace'; %out of the plot list
AVG.Params.defaultPlotName='ActivityTrace'; %out of the plot list
AVG.Params.scrollBarPixels=20;
AVG.Params.shortTimeFrameMovement=0.2; % [0-1] movement relative to window size
AVG.Params.AVversion='1.0';
AVG.Params.scroolYZoomFactor=0.02; %[0-1]
AVG.Params.mouseBottonPressIndicator=0; %an indicator specifying if the mouse botton is currently pressed
AVG.Params.defaultStartTime=0; %[ms]
AVG.Params.defaultWindow=500; %[ms]
test=[]; %a dummy paramter for testing
AVG.Params.exportFigNumber=1; %initialize export figure counter
AVG.Params.maxNumberOfTriggerInColumn=5; %the maximal number of trigger checkbox objects in a column - above this number stacking is also done over rows 
AVG.Params.defaultTriggerOffset=0; % the default value for trigger offset
AVG.Params.selectedAnalysis=[]; %initialize analysis options
AVG.Params.loadTriggerDefault=false; %how to inialize the load trigger mode

%filter
AVG.Params.filterDesigns={'ellip','butter','cheby1','cheby2','equiripple','kaiserwin'}; %should be listed in filter object
AVG.Params.defaultFilterValue=1;
AVG.Params.defaultFilterOrder=8;
AVG.Params.defaultPadding=20; %[ms]
AVG.Params.defaultFilterLowPassStopFreq=2100; %[Hz]
AVG.Params.defaultFilterHighPassPassFreq=200; %[Hz]
AVG.Params.defaultFilterLowPassPassFreq=2000; %[Hz]
AVG.Params.defaultFilterHighPassStopFreq=190; %[Hz]
AVG.Params.defaultRippleInPassband=0.5; %Ap
AVG.Params.defaultAttenuationInHighpass=4; %Ast
AVG.Params.defaultAttenuationInLowpass=4; %Ast
AVG.Params.defaultFilterSamplingRatio=1;
AVG.Params.defaultFilterMethodValue=1;
AVG.Params.samplingFrequency(1)=20000;

if nargin~=floor(nargin/2)*2
    AVG.recordingObj=varargin{1};
    varargin=varargin(2:end);
end
%get input data into AVG.Params. ('property','value' format);
for v=1:2:length(varargin)
    eval(['AVG.Params.' varargin{v} '=' 'varargin{v+1};'])
end

createIcons(); %create icons for GUI buttons
%% %%%%%%%%%%%%%%%%%%%%%%%%%% Initialization  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize params
AVG.Params.filterPadding=AVG.Params.defaultPadding;

%construct plot object - this object contains the data and the different plot types are its methods
AVG.plotData=AVPlotData;
%initiate with arbitrary data
AVG.plotData.M=randn(10,1,100);
AVG.plotData.T=1:100;

%collect all plot types
[AVG.Params.PlotMethods,AVG.Params.PlotList]=AVG.plotData.getPlotNames;
AVG.Params.currentPlotName=AVG.Params.defaultPlotName;

%collect all analysis methods
AVG.Params.analysisMethods=dir([AVG.Params.dataAnalysisDirectory 'analysis*.m']);
AVG.Params.analysisMethods={AVG.Params.analysisMethods.name};
AVG.Params.analysisMethods=cellfun(@(x) x(1:end-2),AVG.Params.analysisMethods,'UniformOutput',0);
AVG.Params.analysisMethodsNames=cellfun(@(x) x(9:end),AVG.Params.analysisMethods,'UniformOutput',0);
%reorder methods so that Filter is the first method
pFilter=strcmp(AVG.Params.analysisMethods,'analysisFilter');
AVG.Params.analysisMethods=[AVG.Params.analysisMethods(pFilter) AVG.Params.analysisMethods(~pFilter)];
AVG.Params.analysisMethodsNames=[AVG.Params.analysisMethodsNames(pFilter) AVG.Params.analysisMethodsNames(~pFilter)];

%get all dataRecording objects
AVG.Params.dataRecordingClasses=dir([AVG.Params.dataClassDirectory '@*']);
AVG.Params.dataRecordingClasses={AVG.Params.dataRecordingClasses.name};
AVG.Params.dataRecordingClasses=cellfun(@(x) x(2:end),AVG.Params.dataRecordingClasses,'UniformOutput',0);

%initialize viewer
AVG.Params.nTriggers=0;
AVG.Params.trigger={};

%construct filter object
AVG.filterObj=filterData; %filter constructor - initialized without sampling frequency
[AVG.Params.filterMethods,AVG.Params.filterNames]=AVG.filterObj.getFilters; %get all filtering methods

%create GUI
createAVGUI;

set(AVG.hMainFigure.hMainAxis,'YLimMode','manual');

if ~isempty(AVG.recordingObj) %initiate viewer with input dataRecording object
    initializeNewRecording;
end

%if timeWindowExternalHandles are passed as input, GUI acts as a slave that is updated everytime the start time or window of the master
%GUI is updated. The color of slave GUI is changed to
if ~isempty(AVG.Params.timeWindowExternalHandles)
    AVG.hSyncStartTime=addlistener(AVG.Params.timeWindowExternalHandles(1),'String','PostSet',@CallbackSyncGUI);
    AVG.hSyncWindow=addlistener(AVG.Params.timeWindowExternalHandles(2),'String','PostSet',@CallbackSyncGUI);
    AVG.hNav.WindowTxt.BackgroundColor=[0.7 1 0.7];AVG.hNav.TimeTxt.BackgroundColor=[0.7 1 0.7];
end
        %% %%%%%%%%%%%%%%%% Nested functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% Update plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function updatePlot()
        %status==0 basic update
        %status==1 channel update
        set(AVG.hGen.messageBox,'string','Plotting data','ForegroundColor','r');drawnow;
        %set channel data
        AVG.plotData.channelNumbers=AVG.Params.channelNumbers(AVG.Params.activeChannelPlaces);
        AVG.plotData.channelNames=AVG.Params.channelNames(AVG.Params.activeChannelPlaces);
        %set activity data
        if any(AVG.Params.selectedAnalysis==1) % a filter was selected -> data is padded
            [AVG.plotData.M,AVG.plotData.T]=AVG.recordingObj.getData(AVG.Params.channelNumbers(AVG.Params.activeChannelPlaces),...
                AVG.Params.startTime-AVG.Params.filterPadding,AVG.Params.window+2*AVG.Params.filterPadding);
            for i=AVG.Params.selectedAnalysis
                if i==1 %case filter
                    [AVG.plotData.M,AVG.plotData.T]=AVG.filterObj.getFilteredData(AVG.plotData.M);
                else
                    feval(AVG.Params.analysisMethods{i},AVG.plotData);
                end
            end
            %implement a cut edges method for the parent class
            AVG.plotData.M=AVG.plotData.M(:,:,round(AVG.Params.paddingSamples/AVG.filterObj.downSamplingFactor) +1:end-round(AVG.Params.paddingSamples/AVG.filterObj.downSamplingFactor));
            AVG.plotData.T=AVG.plotData.T(1:end-round(AVG.Params.paddingSamples/AVG.filterObj.downSamplingFactor)*2); %the first time is 0
        else
            [AVG.plotData.M,AVG.plotData.T]=AVG.recordingObj.getData(...
                AVG.Params.channelNumbers(AVG.Params.activeChannelPlaces),AVG.Params.startTime,AVG.Params.window);
            for i=AVG.Params.selectedAnalysis
                if i==1 %case filter
                    AVG.plotData.M=AVG.filterObj.getFilteredData(AVG.plotData.M);
                else
                    feval(AVG.Params.analysisMethods{i},AVG.plotData);
                end
            end
        end

        if AVG.plotData.plotAnalogChannels
            AVG.plotData.A=AVG.recordingObj.getAnalogData(AVG.Params.analogChannelNumbers(AVG.Params.activeAnalogPlaces),AVG.Params.startTime,AVG.Params.window);
            AVG.plotData.analogChannelNames=AVG.Params.analogChannelNames(AVG.Params.activeAnalogPlaces);
        end

        if AVG.hTrigger.hPlotTrigs.Value && AVG.Params.nCurrentTriggers>0
            AVG.plotData.trigMarks=AVG.Params.currentTriggerSorted(AVG.Params.currentTriggerSorted>AVG.Params.startTime(1) & AVG.Params.currentTriggerSorted<=(AVG.Params.startTime(1)+AVG.Params.window))-AVG.Params.startTime(1);
        end

        %Run the relevant plot method
        AVG.plotData.replot;

        %set slider
        set(AVG.hMainFigure.hTimeSlider,'value',(AVG.Params.startTime(1)+AVG.Params.window) /(AVG.recordingObj.recordingDuration_ms+AVG.Params.window));
        
        set(AVG.hGen.messageBox,'string','Ready','ForegroundColor','k');
    end
%% %%%%%%%%%%%%%%%%%%%%%%%%%% Initialize viewer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function initializeViewer()
        %Recording object
        if iscell(AVG.recordingObj.recordingName)
            AVG.recordingObj.recordingName=AVG.recordingObj.recordingName{1};
        end
        if iscell(AVG.recordingObj.recordingDir)
            AVG.recordingObj.recordingDir=AVG.recordingObj.recordingDir{1};
        end
        if isempty(AVG.recordingObj.recordingName)
            AVG.recordingObj.recordingName='';
        end
        set(AVG.hMainFigure.hFigure,'name',['Activity viewer - ' AVG.recordingObj.recordingDir ' - ' AVG.recordingObj.recordingName]);
        AVG.Params.startDate=AVG.recordingObj.startDate; %show with: datestr(AVG.Params.startDate)
        AVG.Params.endDate=AVG.recordingObj.endDate;
        AVG.Params.recordingTimeLimits=[0 AVG.recordingObj.recordingDuration_ms];
        
        initializeTriggers;
        
        %adjust sampling frequency related features
        AVG.Params.samplingFrequency=AVG.recordingObj.samplingFrequency;
        AVG.Params.Bin=1/AVG.Params.samplingFrequency(1)*1000; %convert from sampling freq [Hz] to time [ms]
        
        %adjust channels
        AVG.Params.channelNames=AVG.recordingObj.channelNames;
        AVG.Params.channelNumbers=AVG.recordingObj.channelNumbers;
        AVG.Params.analogChannelNames=AVG.recordingObj.analogChannelNames;
        AVG.Params.analogChannelNumbers=AVG.recordingObj.analogChannelNumbers;
        
        AVG.Params.activeChannelPlaces=1:numel(AVG.Params.channelNumbers);
        AVG.Params.activeAnalogPlaces=true(1,numel(AVG.Params.analogChannelNumbers));
        
        %adjust start time
        AVG.Params.startTime=AVG.Params.defaultStartTime; %[ms]
        AVG.Params.window=AVG.Params.defaultWindow; %[ms]
        AVG.Params.multiStartTime=numel(AVG.Params.startTime)>1;
        set(AVG.hNav.startTimeEdit,'string',AVG.Params.defaultStartTime);
        set(AVG.hNav.WindowEdit,'string',AVG.Params.defaultWindow);
        set(AVG.hNav.TimeTxt,'string',['start time [ms] /' num2str(AVG.Params.recordingTimeLimits(2))])
        
        %construct filter object
        %AVG.filterObj=filterData(AVG.Params.samplingFrequency(1)); %filter constructor
        AVG.filterObj.samplingFrequency=AVG.Params.samplingFrequency(1);
        set(AVG.hFilter.hSamplingFrequencyEdit,'string',num2str(AVG.Params.samplingFrequency(1)));
        AVG.Params.paddingSamples=AVG.Params.filterPadding/1000*AVG.Params.samplingFrequency(1);
        designFilter; %design the filter
                
        %delete previous plot
        AVG.plotData.deletePlotBackground;
        AVG.plotData.deletePlotControls;
        AVG.plotData.recordingObjReference=AVG.recordingObj; %in case something should be changed in the recording object from the GUI
        
        %initialize new plot
        AVG.plotData.initializePlot(AVG.Params.currentPlotName,AVG.hMainFigure.hMainAxis,AVG.hPlotProp.hMainPanel);
        set(AVG.hGen.messageBox,'string','Ready','ForegroundColor','k');
    end

    function initializeTriggers()
        %adjust trigger related features
        set(AVG.hGen.messageBox,'string','Loading triggers','ForegroundColor','r');drawnow;
        
        if AVG.hTrigger.hGetTrigFromRawFiles.Value
            AVG.Params.triggers=AVG.recordingObj.getTrigger(); %this is the most time consuming step
        else
            AVG.Params.triggers={};
        end
        
        isTriggerActive=cellfun(@(x) ~isempty(x), AVG.Params.triggers); %empty triggers are automatically set to non-active
        if sum(isTriggerActive)>0
            AVG.Params.triggers=AVG.Params.triggers(isTriggerActive);
            AVG.Params.nTriggers=sum(isTriggerActive);
            AVG.Params.currentTrigger=AVG.Params.triggers{1}; %set the first trigger as the current trigger
            [AVG.Params.currentTriggerSorted,AVG.Params.currentTriggerSortedOrder]=sort(AVG.Params.currentTrigger);
            AVG.Params.nCurrentTriggers=numel(AVG.Params.currentTrigger);
        else
            AVG.Params.nTriggers=0;
            AVG.Params.currentTrigger=[];
            AVG.Params.nCurrentTriggers=0;
            AVG.Params.triggers={};
        end
        AVG.Params.triggerOffset=AVG.Params.defaultTriggerOffset;
        createTriggerGUI(); %also updates gui in the case of no triggers
        set(AVG.hGen.messageBox,'string','Ready','ForegroundColor','k');
    end

    function success=downloadGUILayout()
        success=false;
        disp('GUI Layout add on should be installed first!, go to');
        reply = input('Do you want to download it now? Y/N [Y]:','s');
        if isempty(reply)
            reply = 'N';
        end
        if ~strcmpi(reply, 'Y')
            disp('Please install the package and try again.')
            return;
        else
            
            remoteZipFile = 'http://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/47982/versions/7/download/zip';
            try
                saveIn=[AVG.Params.timeSeriesViewerMainDir filesep 'GUILayout'];
                mkdir(saveIn);
                disp(['Downloading toolbox from: ' remoteZipFile]);
                unzip(remoteZipFile, saveIn);
                disp('Download successful! Adding toolbox to the path...');
                % if not in matlab path, add to path
                addpath(genpath(saveIn));
                savepath;
                success=true;
            catch
                error('ERROR: An unexpected error occurred while downloading GUI layout, please install manually')
            end
        end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% General GUI Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function closeMainGUIFigure(hObj,event)
        if ishandle(AVG.hVideoSyncFigure.hFigure)
            close(AVG.hVideoSyncFigure.hFigure); %using close evokes the close request function and is better than delete()
        end
        delete(AVG.plotData);
        delete(AVG.hMainFigure.hFigure);
        clear AVG;
    end
    function CallbackTimeSlider(hObj,event)
        sliderVal=get(AVG.hMainFigure.hTimeSlider,'value');
        AVG.Params.startTime=sliderVal*(AVG.recordingObj.recordingDuration_ms+AVG.Params.window)-AVG.Params.window;
        set(AVG.hNav.startTimeEdit,'string',num2str(AVG.Params.startTime));
        updatePlot;
    end
    function CallbackRestartPlotPush(hObj,event)
        delete(allchild(AVG.hMainFigure.hFigure));
        AVG.plotData.delete;
        AVG.plotData=AVPlotData;
        createAVGUI(AVG.hMainFigure.hFigure);
        createTriggerGUI();
        
        set(AVG.hMainFigure.hMainAxis,'YLimMode','manual');
        
        initializeViewer;
        updatePlot;
    end
    function CallbackRefreshPlotPush(hObj,Event)
        updatePlot;
    end
    function CallbackExportPlotPush(hObj,Event)
        cmap=get(AVG.hMainFigure.hMainAxis,'colormap');
        hNewFig=figure('Position',[10 65 1700 940]);
        copyobj(AVG.hMainFigure.hMainAxis,hNewFig);
        
        hAxes=get(hNewFig,'Children');
        hAxes(1).Units='normalized';
        set(hAxes(1),'FontSize',12,'Position',[0.07 0.07 0.91 0.91]);
        set(get(hAxes(1),'xlabel'),'fontsize',14);
        set(get(hAxes(1),'ylabel'),'fontsize',14);
        
        set(hNewFig,'Colormap',cmap,'PaperPositionMode','auto');
        assignin('base',['hExportedPlot' num2str(AVG.Params.exportFigNumber)],hNewFig);
        AVG.Params.exportFigNumber=AVG.Params.exportFigNumber+1;
    end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Menu Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function CallbackExportData(hObj,event)
        dlgOptions.WindowStyle='normal';
        dataVarible = inputdlg({'Please enter variable name:'},'Export data to workspace',1,{''},dlgOptions);
        assignin('base',[dataVarible{1} '_T'],AVG.plotData.T);
        assignin('base',[dataVarible{1} '_M'],AVG.plotData.M);
        disp('In exported data, M is structured as: [Channels x Trials x Time] and T is the corresponding times in [ms]');
    end

    function addSyncAV(hObj,event)
        timeSeriesViewer('timeWindowExternalHandles',[AVG.hNav.startTimeEdit AVG.hNav.WindowEdit]);
    end

    function addSyncVideo(hObj,event)
        %set default params
        AVG.Params.videoSyncVerified=false;
        AVG.Params.currentSpeedup=1;

        if isempty(AVG.hVideoSyncFigure.hFigure) || ~isvalid(AVG.hVideoSyncFigure.hFigure)
            createSyncVideoGUI;
        else
            figure(AVG.hVideoSyncFigure.hFigure)
        end
    end

    function closeSyncedVideoFigure(hObj,event)
        delete(AVG.hVideoSyncFigure.hFigure);
        delete(AVG.Params.videoReader);
        %close video player objects
    end

    function CallbackLoadVideoEdit(hObj,event)
        AVG.hVideoSyncFigure.hCheckSyncPush.BackgroundColor=[0.8 0.8 0.8]; %bring background to origin to show that sync was not checked
        AVG.Params.videoSyncVerified=false;
    end

    function CallbackLoadVideoPush(hObj,event)
        if isempty(AVG.recordingObj)
            startDir=cd;
        else
            startDir=AVG.recordingObj.recordingDir;
        end
        [tmpVideoFile,tmpVideoDir]= uigetfile('*.*','Choose the video file',startDir);
        if tmpVideoFile==0 %no folder chosen
            msgbox('File was not chosen. Video were not added!','Attention','error','replace');
        else
            AVG.hVideoSyncFigure.hLoadVideoEdit.String=fullfile(tmpVideoDir,tmpVideoFile);
        end
        AVG.hVideoSyncFigure.hCheckSyncPush.BackgroundColor=[0.8 0.8 0.8]; %bring background to origin to show that sync was not checked
        AVG.Params.videoSyncVerified=false;
    end

    function CallbackCheckSyncPush(hObj,event)
        if AVG.Params.nTriggers==0
            msgbox('No triggers exist! Please load trigger and try again','Attention','error','replace');
        else
            %!!!! Check if can read grey scale video instead of RGB to save space and time ('Grayscale')
            AVG.Params.videoReader   = VideoReader(AVG.hVideoSyncFigure.hLoadVideoEdit.String);
            AVG.Params.videoDuration=AVG.Params.videoReader.Duration;
            AVG.Params.frameRate=AVG.Params.videoReader.FrameRate;
            AVG.Params.frameRateOriginal=AVG.Params.frameRate; % in some cases the video meta data is not correct. For some calculation the original value is still be important.
            AVG.Params.frameRate=1000/median(diff(AVG.Params.currentTrigger(1:100)));
            AVG.hVideoSyncFigure.hFrameRateEdit.String=num2str(AVG.Params.frameRate);
            AVG.Params.nFramesVideo=round(AVG.Params.videoDuration*AVG.Params.frameRateOriginal);

            AVG.Params.triggerFrameSync=AVG.Params.currentTrigger;

            %checking for 1 sec trigger gaps indicating start and stop recording session
            nFrameGaps=numel(find(diff(AVG.Params.currentTrigger)>=1000));
            if nFrameGaps>=2
                if nFrameGaps==2 %this may happen when recording with reptiLearn
                    fprintf('Exactly two >1 sec gaps found in triggers. Removing irrelevant triggers before and after gaps\n');
                elseif nFrameGaps>2
                    fprintf('More than two >1 sec gaps found in triggers. Assuming relevant frames are within the largest gap...\n');
                end
                pStartEnd=find(diff(AVG.Params.triggerFrameSync)>1000 | diff(AVG.Params.triggerFrameSync)<-1000);
                [~,maxDiff]=max(diff(pStartEnd));
                pStartEnd=pStartEnd(maxDiff:maxDiff+1);
                AVG.Params.triggerFrameSync=AVG.Params.triggerFrameSync((pStartEnd(1)+1):pStartEnd(2));
            elseif nFrameGaps<=1
                fprintf('One or zero >1 sec gaps found in triggers... Using all triggers for analysis!\n');
            end

            AVG.Params.diffFrames=numel(AVG.Params.triggerFrameSync)-round(AVG.Params.nFramesVideo);
            AVG.Params.pSync=[];% the relevant places in AVG.Params.triggerFrameSync
            if AVG.Params.diffFrames==0
                AVG.hVideoSyncFigure.hCheckSyncPush.BackgroundColor=[0 1 0];
            elseif AVG.Params.diffFrames>0
                if all(AVG.hVideoSyncFigure.hValidSyncTriggersEdit.String)=='0'
                    msgbox({['Found ' num2str(AVG.Params.diffFrames) ' more frames in digital triggers than in video (' num2str(round(AVG.Params.nFramesVideo)) ')!!!'],...
                        'Proceeding with analysis assuming uniform distribution of lost frames in video'},'Attention','error','replace');
                    p2Remove=true(1,numel(AVG.Params.triggerFrameSync));
                    p2Remove(round((1:AVG.Params.diffFrames)/AVG.Params.diffFrames*numel(AVG.Params.triggerFrameSync)))=false;
                    AVG.Params.pSync=find(p2Remove);
                else
                    msgbox({['Found ' num2str(AVG.Params.diffFrames) ' more frames in digital triggers than in video (' num2str(round(AVG.Params.nFramesVideo)) ')!!!'],...
                        'Proceeding with analysis assuming frames 1:numberOfFrames'},'Attention','error','replace');
                    AVG.hVideoSyncFigure.hValidSyncTriggersEdit.String=['1:' num2str((numel(AVG.Params.triggerFrameSync)-AVG.Params.diffFrames))];
                    AVG.hVideoSyncFigure.hCheckSyncPush.BackgroundColor=[0.5 0 0];
                    AVG.Params.pSync=str2num(AVG.hVideoSyncFigure.hValidSyncTriggersEdit.String);
                end
            end
            AVG.Params.videoSyncVerified=true;
        end
    end

    function CallbackConvertVideo2EphysEdit(hObj,event)
        timeStr=hObj.String;
        [~, ~, D, H, MN, S] = datevec(timeStr);
        frameNum=round(H*60*60+MN*60+S)*AVG.Params.frameRateOriginal;
        ephysTime=AVG.Params.triggerFrameSync(AVG.Params.pSync(frameNum));
        AVG.hVideoSyncFigure.hConvertEphys2VideoEdit.String=num2str(ephysTime);
    end

    function CallbackConvertEphys2VideoEdit(hObj,event)
        tStr=str2num(hObj.String);
        t=find(AVG.Params.triggerFrameSync(AVG.Params.pSync)>tStr,1,'first')/AVG.Params.frameRateOriginal*1000;
        if ~isempty(t)
            h=floor(t/(60*60*1000));
            m=floor((t-h*60*60*1000)/(60*1000));
            s=floor((t-h*60*60*1000-m*60*1000)/1000);
            ms=floor(t-h*60*60*1000-m*60*1000-s*1000);
            timeStr=num2str(h)+":"+num2str(m)+":"+num2str(s)+"."+num2str(ms);
            AVG.hVideoSyncFigure.hConvertVideo2EphysEdit.String=timeStr;
        else
            AVG.hVideoSyncFigure.hConvertVideo2EphysEdit.String='Unvalid time!';
        end
    end

    function CallbackExportIMPlayPush(hObj,event)
        hObj.BackgroundColor=[1 0 0];

        if AVG.Params.videoSyncVerified
            pFrames=find(AVG.Params.triggerFrameSync(AVG.Params.pSync)>AVG.Params.startTime & AVG.Params.triggerFrameSync(AVG.Params.pSync)<(AVG.Params.startTime+AVG.Params.window));
            nActFrames=numel(pFrames);
            nFrames=AVG.Params.window/1000*AVG.Params.frameRate;
            if nActFrames<nFrames-1 %can happen if the there is no video (of full video during the relevant times
                msgbox(['Missing frames in segment! video exists between ' num2str(AVG.Params.triggerFrameSync(AVG.Params.pSync(1))) ' - ' num2str(AVG.Params.triggerFrameSync(AVG.Params.pSync(end)))],'Attention','error','replace');
                return;
            end
            AVG.Params.videoReader.CurrentTime = pFrames(1)/AVG.Params.frameRateOriginal;
            frameVid=AVG.Params.videoReader.readFrame;
            videoMat=zeros([size(frameVid) nActFrames],class(frameVid));
            videoMat(:,:,:,1)=frameVid;
            for frames=2:nActFrames
                videoMat(:,:,:,frames)=AVG.Params.videoReader.readFrame; %generalize to grey scale video.
            end
            implay(videoMat,AVG.Params.frameRate);
        else
            msgbox('Triggers not synced to video, run sync first','Attention','error','replace');
        end
        hObj.BackgroundColor=[0.8 0.8 0.8];
    end

    function CallbackPlaySpeedupPush(hObj,event)
        AVG.Params.currentSpeedup=str2num(hObj.String(2:end));
        AVG.Params.currentSpeedup=4*AVG.Params.currentSpeedup;
        if AVG.Params.currentSpeedup>32
            AVG.Params.currentSpeedup=1;
        end
        hObj.String=['X' num2str(AVG.Params.currentSpeedup)];
    end

    function CallbackRunVideoPush(hObj,event)
        hObj.BackgroundColor=[0 1 0];
        if AVG.Params.videoSyncVerified
            pFrames=find(AVG.Params.triggerFrameSync(AVG.Params.pSync)>AVG.Params.startTime & AVG.Params.triggerFrameSync(AVG.Params.pSync)<(AVG.Params.startTime+AVG.Params.window));
            %startFrame=find(AVG.Params.triggerFrameSync(AVG.Params.pSync)>AVG.Params.startTime,1,'first');
            %endFrame=find(AVG.Params.triggerFrameSync(AVG.Params.pSync)<(AVG.Params.startTime+AVG.Params.window),1,'last');
            nActFrames=numel(pFrames);
            nFrames=AVG.Params.window/1000*AVG.Params.frameRate;
            if nActFrames<nFrames-1 %can happen if the there is no video (of full video during the relevant times
                fprintf('Warning: Missing %d frame triggers in segment or sampling rate missmatch!\nNotice that video exists between %fms-%fms\n',nFrames-nActFrames,AVG.Params.triggerFrameSync(AVG.Params.pSync(1)),AVG.Params.triggerFrameSync(AVG.Params.pSync(end)));
            end
            %% 
            hTmp=line(AVG.hMainFigure.hMainAxis,[0 0],AVG.hMainFigure.hMainAxis.YLim,'color','r');
            hold(AVG.hVideoSyncFigure.hVideoAxis,'on');
            set(AVG.hVideoSyncFigure.hVideoAxis,'nextplot','replacechildren');
            AVG.Params.videoReader.CurrentTime = pFrames(1)/AVG.Params.frameRateOriginal;
            for frames=1:AVG.Params.currentSpeedup:nActFrames
                frameVid=AVG.Params.videoReader.readFrame; %generalize to grey scale video.
                %image(frameVid,'Parent',AVG.hVideoSyncFigure.hVideoAxis);
                imshow(frameVid,'Parent',AVG.hVideoSyncFigure.hVideoAxis);
                tmpTime=AVG.Params.triggerFrameSync(AVG.Params.pSync(pFrames(frames)))-AVG.Params.startTime;
                hTmp.XData=[tmpTime tmpTime];
                drawnow nocallbacks;%nocallbacks;limitrate
                for skipFrame=1:(AVG.Params.currentSpeedup-1)
                    AVG.Params.videoReader.readFrame;
                end
            end
            delete(hTmp);
            hold(AVG.hVideoSyncFigure.hVideoAxis,'off');
        else
            msgbox('Triggers not synced to video, run sync first','Attention','error','replace');
        end
        hObj.BackgroundColor=[0.8 0.8 0.8];
    end

    function CallbackValidSyncTriggersEdit(hObj,event)
        AVG.Params.pSync=str2num(AVG.hVideoSyncFigure.hValidSyncTriggersEdit.String);
        if AVG.hVideoSyncFigure.hValidSyncTriggersEdit.String(1)=='0' %for uniformly distributing triggers
            CallbackCheckSyncPush;
        end
    end

    function CallbackFrameRateEdit(hObj,event)
        AVG.Params.frameRate=str2num(hObj.String);
    end



    function CallbackSyncGUI(hObj,event)
        AVG.Params.window=str2double(AVG.Params.timeWindowExternalHandles(2).String);
        AVG.hNav.startTimeEdit.String=AVG.Params.timeWindowExternalHandles(1).String;
        CallbackNavTimeEdit();
    end

    function CallbackDataTypeSelection(hObj,event)
        %set all fields in menu to off
        for i=1:length(AVG.Params.dataRecordingClasses)
            set(AVG.hMenubar.hRecording.(AVG.Params.dataRecordingClasses{i}),'checked','off');
        end
        set(hObj,'Checked','on');
        AVG.Params.currentDataType=get(hObj,'Label');
        
        eval(['AVG.recordingObj=' AVG.Params.currentDataType ';']); %create a data recording object
        initializeNewRecording;
    end

    function CallbackCredits(hObj,event)
        msgbox({'Developed by Mark Shein-Idelson','e-mail: shein.mark@gmail.com'},'Main developers and contributers');
    end

    function initializeNewRecording
        set(AVG.hGen.messageBox, 'string','Initializing new data','ForegroundColor','r');
        initializeViewer;
        if ishandle(AVG.hVideoSyncFigure.hFigure)
            close(AVG.hVideoSyncFigure.hFigure); %using close evokes the close request function and is better than delete()
        end
        AVG.plotData.refreshPlot=1;
        updatePlot;
    end

    function CallbackPlotTypeSelection(hObj,event)
        %set all fields in menu to off
        for i=1:numel(AVG.Params.PlotMethods)
            set(AVG.hMenubar.hPlot.(AVG.Params.PlotList{i}),'checked','off');
        end
        set(hObj,'Checked','on');
        AVG.Params.currentPlotName=get(hObj,'Label');
        
        %delete previous plot
        AVG.plotData.deletePlotBackground;
        AVG.plotData.deletePlotControls;
        
        AVG.plotData.initializePlot(AVG.Params.currentPlotName,AVG.hMainFigure.hMainAxis,AVG.hPlotProp.hMainPanel);
        AVG.plotData.refreshPlot=1;
        updatePlot;
    end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mouse Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for the mouse callback to be activated the text below should be entered into section creating the GUI
%        iptPointerManager(AVG.hMainFigure.hFigure);
%        mainPlotPointerBehaviour.enterFcn=@CallbackMouseEnterMainPlot;
%        mainPlotPointerBehaviour.exitFcn=@CallbackMouseExitMainPlot;
%        mainPlotPointerBehaviour.traverseFcn=@CallbackMouseTransverseMainPlot;
%        iptSetPointerBehavior(AVG.hMainFigure.hMainAxis, mainPlotPointerBehaviour);
%in addition, the following properties should be set in the figure handle:
%        'WindowButtonMotionFcn',@CallbackMouseMotion,'WindowButtonDownFcn',@CallbackMouseBottonPress, 'WindowButtonUpFcn',@CallbackMouseBottonRelease,'WindowScrollWheelFcn',@CallbackMouseWheelScroll


    function CallbackMouseMotion(hObj,event)
        %AVG.mouse.Position=get(hObj, 'currentpoint'); % get mouse location on figure
    end
    function CallbackMouseWheelScroll(hObj,event)
        yl=ylim(AVG.hMainFigure.hMainAxis);
        dist=yl(2)-yl(1);
        if AVG.Params.mouseBottonPressIndicator==0
            ylim(AVG.hMainFigure.hMainAxis,yl+[-1 1]*dist*AVG.Params.scroolYZoomFactor*event.VerticalScrollCount*event.VerticalScrollAmount)
        else
            ylim(AVG.hMainFigure.hMainAxis,yl+[1 1]*dist*AVG.Params.scroolYZoomFactor*event.VerticalScrollCount*event.VerticalScrollAmount)
        end
    end
    function CallbackMouseBottonPress(hObj,event)
        AVG.Params.mouseBottonPressIndicator=1;
    end
    function CallbackMouseBottonRelease(hObj,event)
        AVG.Params.mouseBottonPressIndicator=0;
    end
    function CallbackMouseEnterMainPlot(hFig,currentPoint)
        set(hFig,'Pointer','hand');
    end
    function CallbackMouseExitMainPlot(hFig,currentPoint)
        set(hFig,'Pointer','arrow');
    end
    function CallbackMouseTransverseMainPlot(hFig,currentPoint)
    end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Navigation Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function CallbackNavTimeFrameMovePush(hObj,event,mov)
        newStartTime=AVG.Params.startTime+mov*AVG.Params.window;
        if any(newStartTime>=(AVG.Params.recordingTimeLimits(2)))
            msgbox('End of recording reached','Attention','error','replace');
        elseif any(newStartTime<AVG.Params.recordingTimeLimits(1)-AVG.Params.window)
            msgbox('Beginning of recording reached','Attention','error','replace');
        else
            AVG.Params.startTime=newStartTime;
            set(AVG.hNav.startTimeEdit,'string',num2str(newStartTime));
            updatePlot;
        end
    end

    function CallbackNavZoomPush(hObj,event,windowChange)
        winCenter=AVG.Params.startTime+AVG.Params.window/2;
        AVG.Params.window=AVG.Params.window*windowChange;
        AVG.Params.startTime=winCenter-AVG.Params.window/2;
        set(AVG.hNav.WindowEdit,'string',num2str(AVG.Params.window));
        set(AVG.hNav.startTimeEdit,'string',num2str(AVG.Params.startTime));
        updatePlot;
    end

    function CallbackNavTimeEdit(hObj,event)
        newStartTime=str2num(get(AVG.hNav.startTimeEdit,'string'));
        if any(newStartTime>AVG.Params.recordingTimeLimits(2))
        %if any(newStartTime>(AVG.Params.recordingTimeLimits(2)-2*AVG.Params.window)) %changed to above on 091014 for testing
            msgbox('End of recording reached','Attention','error','replace');
            set(AVG.hNav.startTimeEdit,'string',num2str(AVG.Params.startTime));
        elseif any(newStartTime<AVG.Params.recordingTimeLimits(1)-AVG.Params.window)
            msgbox('Beginning of recording reached','Attention','error','replace');
            set(AVG.hNav.startTimeEdit,'string',num2str(AVG.Params.startTime));
        else
            if numel(newStartTime)<=1 && AVG.Params.multiStartTime
                set(AVG.hNav.TimeTxt,'string','start time [ms]','ForegroundColor','k');
                AVG.Params.multiStartTime=0;
                AVG.plotData.refreshPlot=1;
            elseif numel(newStartTime)>1
                set(AVG.hNav.TimeTxt,'string','start times [ms] (multi-mode)','ForegroundColor','r');
                AVG.Params.multiStartTime=1;
            end
            AVG.Params.startTime=newStartTime;
            updatePlot;
        end
    end

    function CallbackNavWindowEdit(hObj,event)
        AVG.Params.window=str2num(get(hObj,'string'));
        updatePlot;
    end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Channel selection Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Todo: All the channel selection functions can be simplified for the electrode channel if activeChannelPlaces is a binary vector. See how it is done for the analog channels
    function CallbackChAddChannelsEdit(hObj,event)
        inStr=get(hObj,'string');
        ChannelsToAdd=str2num(inStr);
        if ~isempty(ChannelsToAdd)
            for i=1:numel(ChannelsToAdd)
                P=find(AVG.Params.channelNumbers==ChannelsToAdd(i));
                if ~isempty(P)
                    AVG.Params.activeChannelPlaces=[AVG.Params.activeChannelPlaces P];
                end
            end
        else
            AVG.Params.activeAnalogPlaces=AVG.Params.activeAnalogPlaces | ismember(AVG.Params.analogChannelNames,split(inStr,' '));
        end
        set(hObj,'string','Add channels'); %replace entered channel string with default string
        updatePlot;
    end

    function CallbackChRemoveChannelsEdit(hObj,event)
        %If a channel that was selected does not exist or is already removed, it is not removed and no error is given
        inStr=get(hObj,'string');
        ChannelsToRemove=str2num(inStr);
        if ~isempty(ChannelsToRemove)
            for i=1:length(ChannelsToRemove)
                AVG.Params.activeChannelPlaces(find(AVG.Params.channelNumbers(AVG.Params.activeChannelPlaces)==ChannelsToRemove(i)))=[];
            end
        else
            AVG.Params.activeAnalogPlaces(ismember(AVG.Params.analogChannelNames,split(inStr,' ')))=false;
        end
        set(hObj,'string','Remove channels');%replace entered channel string with default string
        updatePlot;
    end

    function CallbackChKeepOnlyChannelsEdit(hObj,event)
        inStr=get(hObj,'string');
        ChannelsToKeep=str2num(inStr);
        if ~isempty(ChannelsToKeep)
            tmpActiveChannels=zeros(1,numel(ChannelsToKeep));
            for i=1:length(ChannelsToKeep)
                P=find(AVG.Params.channelNumbers==ChannelsToKeep(i));
                if ~isempty(P)
                    tmpActiveChannels(i)=P;
                else
                    disp('Part of the selected channels do not exist in the recording!!!!');
                    return;
                end
            end
            AVG.Params.activeChannelPlaces=tmpActiveChannels;
        else
            AVG.Params.activeAnalogPlaces=ismember(AVG.Params.analogChannelNames,split(inStr,' '));
        end
        set(hObj,'string','Keep only channels');%replace entered channel string with default string
        updatePlot;
    end

    function CallbackChEnableAllChannelsPush(hObj,event)
        AVG.Params.activeChannelPlaces=1:numel(AVG.Params.channelNumbers);
        AVG.Params.activeAnalogPlaces=true(1,numel(AVG.Params.analogChannelNumbers));
        updatePlot;
    end

    function CallbackChManualChannelSelectPush(hObj,event)
        channelsPerColumn=30;
        nActiveChannels=length(AVG.Params.activeChannelPlaces);
        nChannels=numel(AVG.Params.channelNumbers);
        nActiveAnalogChannels=length(AVG.Params.activeAnalogPlaces);
        nAnalogChannels=numel(AVG.Params.analogChannelNumbers);
        
        AVG.manualChannelSelectionFigure.hFigure = figure('Position',[AVG.hMainFigure.scrsz(3)*0.01 AVG.hMainFigure.scrsz(4)*0.07 AVG.hMainFigure.scrsz(3)*0.3 AVG.hMainFigure.scrsz(4)*0.85], ...
            'Name','Activity Viewer - manual channel selection', 'NumberTitle','off', 'MenuBar','none', 'Toolbar','none', 'HandleVisibility','off');
        AVG.manualChannelSelectionFigure.VBox=uix.VBox('Parent',AVG.manualChannelSelectionFigure.hFigure, 'Spacing',4, 'Padding',4);
        AVG.manualChannelSelectionFigure.hChannelGrid = uix.Grid('Parent',AVG.manualChannelSelectionFigure.VBox, 'Spacing',5, 'Padding',10);
        AVG.manualChannelSelectionFigure.applySelectionPush=uicontrol('Parent', AVG.manualChannelSelectionFigure.VBox, 'Callback',{@CallbackManualChSelectionApplyPush,}, 'Style','push', 'String','Apply selection');
        for i=1:nChannels
            AVG.manualChannelSelectionFigure.hCh(i)=uicontrol('Parent', AVG.manualChannelSelectionFigure.hChannelGrid, ...
                'Style','checkbox', 'String', [num2str(AVG.Params.channelNumbers(i)) ' - ' AVG.Params.channelNames{i}],'value',0);
        end
        for i=1:nAnalogChannels
            AVG.manualChannelSelectionFigure.hCh(i+nChannels)=uicontrol('Parent', AVG.manualChannelSelectionFigure.hChannelGrid, ...
                'Style','checkbox', 'String', [num2str(AVG.Params.analogChannelNumbers(i)) ' - ' AVG.Params.analogChannelNames{i}],'value',0);
        end
        for i=1:nActiveChannels
            set(AVG.manualChannelSelectionFigure.hCh(AVG.Params.activeChannelPlaces),'value',1);
        end
        for i=1:nActiveAnalogChannels
            set(AVG.manualChannelSelectionFigure.hCh(nChannels+find(AVG.Params.activeAnalogPlaces)),'value',1);
        end
        set(AVG.manualChannelSelectionFigure.hChannelGrid,'Widths',-1*ones(1,ceil((nActiveChannels+nActiveAnalogChannels)/channelsPerColumn)),'Heights', -1*ones(1,channelsPerColumn));
        set(AVG.manualChannelSelectionFigure.VBox, 'Heights',[-1 30]);
    end

    function CallbackManualChSelectionApplyPush(hObj,event)
        nCh=numel(AVG.Params.channelNumbers);
        for i=1:(nCh+numel(AVG.Params.analogChannelNumbers))
            p(i)=get(AVG.manualChannelSelectionFigure.hCh(i),'value');
        end
        AVG.Params.activeChannelPlaces=AVG.Params.channelNumbers(p(1:nCh)==1);
        AVG.Params.activeAnalogPlaces=AVG.Params.analogChannelNumbers(p((nCh+1):end)==1);
        delete(AVG.manualChannelSelectionFigure.hFigure);
        updatePlot;
    end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filter Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function CallbackFilterLowPassPassEdit(hObj,event),designFilter,end
    function CallbackFilterHighPassPassEdit(hObj,event),designFilter,end
    function CallbackFilterLowPassStopEdit(hObj,event),designFilter,end
    function CallbackFilterHighPassStopEdit(hObj,event),designFilter,end
    function CallbackFilterMethodPopup(hObj,event),designFilter,end
    function CallbackFilterOrderEdit(hObj,event),designFilter,end
    function CallbackFilterSamplingRatioEdit(hObj,event),designFilter,end
    function CallbackFilterSamplingFrequencyEdit(hObj,event),designFilter,end
    function CallbackFilterDesignPopup(hObj,event),designFilter,end
    function CallbackFilterRippleEdit(hObj,event),designFilter,end
    function CallbackFilterLowPassAttenuationEdit(hObj,event),designFilter,end
    function CallbackFilterHighPassAttenuationEdit(hObj,event),designFilter,end

    function CallbackPaddingSamplesEdit(hObj,event)
        AVG.Params.filterPadding=str2num(get(hObj,'string'));
        AVG.Params.paddingSamples=AVG.Params.filterPadding/1000*AVG.Params.samplingFrequency(1);
    end
    function CallbackPlotFilterPush(hObj,event)
        AVG.filterObj.plotFilter;
    end
    function designFilter()
        methodValue=get(AVG.hFilter.hMethodPopup,'value');
        AVG.filterObj.filterOrder=str2num(get(AVG.hFilter.hOrderEdit,'string'));
        AVG.filterObj.filterDesign=AVG.Params.filterDesigns{get(AVG.hFilter.hDesignPopup,'value')};
        AVG.filterObj.highPassPassCutoff=str2num(get(AVG.hFilter.hHighPassPassEdit,'string'));
        AVG.filterObj.highPassStopCutoff=str2num(get(AVG.hFilter.hHighPassStopEdit,'string'));
        AVG.filterObj.lowPassPassCutoff=str2num(get(AVG.hFilter.hLowPassPassEdit,'string'));
        AVG.filterObj.lowPassStopCutoff=str2num(get(AVG.hFilter.hLowPassStopEdit,'string'));
        AVG.filterObj.rippleInPassband=str2num(get(AVG.hFilter.hRippleEdit,'string'));
        AVG.filterObj.attenuationInHighpass=str2num(get(AVG.hFilter.hHighpassAttenuationEdit,'string'));
        AVG.filterObj.attenuationInLowpass=str2num(get(AVG.hFilter.hLowpassAttenuationEdit,'string'));
        AVG.filterObj.downSamplingFactor=str2num(get(AVG.hFilter.hSamplingRatioEdit,'string'));
        
        if ~strcmp(AVG.filterObj.currentDesignedFilter,'bandstop')
            if AVG.Params.samplingFrequency(1)<2*AVG.filterObj.lowPassPassCutoff
                AVG.filterObj.lowPassStopCutoff=AVG.Params.samplingFrequency(1)/2;
                AVG.filterObj.lowPassPassCutoff=AVG.filterObj.lowPassStopCutoff*0.95;
                set(AVG.hFilter.hHighPassPassEdit,'string',num2str(AVG.filterObj.highPassPassCutoff));
                set(AVG.hFilter.hHighPassStopEdit,'string',num2str(AVG.filterObj.highPassStopCutoff));
                disp('Notice: filter settings were changed since sampling frequency was lower than the highpass cutoff');
            end
            if AVG.filterObj.lowPassPassCutoff<=AVG.filterObj.highPassPassCutoff
                
                AVG.filterObj.highPassStopCutoff=AVG.filterObj.lowPassPassCutoff/2;
                AVG.filterObj.highPassPassCutoff=AVG.filterObj.highPassStopCutoff*1.05;
                
                set(AVG.hFilter.hHighPassPassEdit,'string',num2str(AVG.filterObj.highPassPassCutoff));
                set(AVG.hFilter.hHighPassStopEdit,'string',num2str(AVG.filterObj.highPassStopCutoff));
                disp('Notice: filter settings were changed since highpass cutoff was higher than lowpass cutoff');
            end
        end
        
        AVG.filterObj=AVG.filterObj.(AVG.Params.filterMethods{methodValue}); %design the filter
        set(AVG.hFilter.hOrderEdit,'string',num2str(AVG.filterObj.filterOrder));
    end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analysis Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function CallbackAnaysisTogglePush(hObj,Event,select)
        if get(hObj,'value')
            AVG.Params.selectedAnalysis=[AVG.Params.selectedAnalysis select];
        else
            set(AVG.hAnalysis.hToggleButtons(select),'string',[AVG.Params.analysisMethodsNames{select} ' (-)']);
            AVG.Params.selectedAnalysis(find(AVG.Params.selectedAnalysis==select))=[];
        end
        for i=1:numel(AVG.Params.selectedAnalysis)
            set(AVG.hAnalysis.hToggleButtons(AVG.Params.selectedAnalysis(i)),'string',[AVG.Params.analysisMethodsNames{AVG.Params.selectedAnalysis(i)} ' (' num2str(i) ')']);
        end
    end
    function CallbackRemoveAllAnalysisPush(hObj,event)
        set(AVG.hAnalysis.hToggleButtons,'value',0);
        for i=1:numel(AVG.hAnalysis.hToggleButtons)
            set(AVG.hAnalysis.hToggleButtons(i),'string',[AVG.Params.analysisMethodsNames{i} ' (-)']);
        end
        AVG.Params.selectedAnalysis=[];
    end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Trigger Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function CallbackManualSetPush(hObj,Event)
        trigIn=get(AVG.hTrigger.manualSetEdit,'string');
        if ~isempty(trigIn)
            tmpTrig={str2num(trigIn)};
            if AVG.Params.nTriggers==0
                AVG.Params.triggers=tmpTrig;
            else
                AVG.Params.triggers=[AVG.Params.triggers tmpTrig];
            end
            AVG.Params.nTriggers=AVG.Params.nTriggers+1;
            createTriggerGUI();
            set(AVG.hTrigger.manualSetEdit,'string','');
            CallbackTriggerCheckbox([],[],AVG.Params.nTriggers)
        else
            disp('Trigger input not provided!');
        end
    end

    function CallbackTrigOffset(hObj,Event)
        tmp=str2num(get(AVG.hTrigger.hOffsetEdit,'string'));
        if isnumeric(tmp) && numel(tmp)==1
            AVG.Params.triggerOffset=tmp;
        else
            msgbox('Input delay is not a valid number','Attention','error','replace');
            set(AVG.hTrigger.hOffsetEdit,'string',0);
        end
    end

    function CallbackManualLoadFile(hObj,Event)
        [tmpTrigFile,tmpTrigDir]= uigetfile('*.mat','Choose the Mat file',cd,'MultiSelect','on');
        if tmpTrigFile==0 %no folder chosen
            msgbox('File was not chosen and triggers were not added','Attention','error','replace');
        else
            tmpTrig=load([tmpTrigDir dirSep tmpTrigFile]);
            trigName=fields(tmpTrig);
            nVars=numel(trigName);
            if nVars==1
                trigName=trigName{1};
            else
                cellVar=[];
                for i=1:nVars
                    if iscell(tmpTrig.(trigName{i}))
                        cellVar=[cellVar i];
                    end
                end
                if numel(cellVar)==1
                    trigName=trigName{cellVar};
                    fprintf('setting triggers to field %s in file since this is the only cell array',trigName);
                elseif numel(cellVar)==0
                    error('Trigger variable should be a cell array of trigger vectors');
                else
                    pTmp=find(strcmp(trigName,{'tTrig'}));
                    if ~isempty(pTmp)
                        trigName=trigName{pTmp};
                    else
                        pTmp=find(strcmp(trigName,{'triggers'}));
                        if ~isempty(pTmp)
                            trigName=trigName{pTmp};
                        else
                            error('No cell arrays named triggers or tTrig were found in file!');
                        end
                    end
                    
                end
            end
            tmpTrig=tmpTrig.(trigName);

            AVG.Params.nTriggers=AVG.Params.nTriggers+numel(tmpTrig);
            AVG.Params.triggers=[AVG.Params.triggers tmpTrig];
            createTriggerGUI();
            CallbackTriggerCheckbox([],[],AVG.Params.nTriggers);
        end
    end

    function CallbackManualLoadVariable(hObj,Event)
        dlgOptions.WindowStyle='normal';
        triggerVar = inputdlg({'Please enter variable name:'},'Get trigger from workspace',1,{''},dlgOptions);
        tmpTrig=evalin('base',triggerVar{1});
        if iscell(tmpTrig)
            AVG.Params.nTriggers=AVG.Params.nTriggers+numel(tmpTrig);
            AVG.Params.triggers=[AVG.Params.triggers tmpTrig];
        elseif isnumeric(tmpTrig)
            AVG.Params.nTriggers=AVG.Params.nTriggers+1;
            AVG.Params.triggers=[AVG.Params.triggers {tmpTrig}];            
        else
            error('Trigger variable should be a cell array of trigger vectors, or a numeric array with one vector');
        end

        createTriggerGUI();
        CallbackTriggerCheckbox([],[],AVG.Params.nTriggers);
    end

    function CallbackManualLoadSpikes(hObj,Event)
        [tmpTrigFile,tmpTrigDir]= uigetfile('*.mat','Choose the Mat file',cd);
        if tmpTrigFile==0 %no folder chosen
            msgbox('File was not chosen and triggers were not added','Attention','error','replace');
        else
            spikeSort=load([tmpTrigDir dirSep tmpTrigFile],'t','ic');
            
            [pSelect] = listdlg('PromptString','Choose a neuron:','SelectionMode','multiple','ListString',num2str(spikeSort.ic(1:2,:)'))
            
            trigName='SS';
            for i=1:numel(pSelect)
                tmpTrig{i}=spikeSort.t(spikeSort.ic(3,pSelect(i)):spikeSort.ic(4,pSelect(i)));
            end

            AVG.Params.nTriggers=AVG.Params.nTriggers+numel(tmpTrig);
            AVG.Params.triggers=[AVG.Params.triggers tmpTrig];
            createTriggerGUI();
            CallbackTriggerCheckbox([],[],AVG.Params.nTriggers);
        end
    end

    function CallbackTriggerNumberEdit(hObj,Event)
        tmpTriggerNumber=str2num(get(AVG.hTrigger.hNumber,'string'));
        if tmpTriggerNumber>AVG.Params.nCurrentTriggers || tmpTriggerNumber<1
            msgbox('Selected trigger number is not within current trigger limits','Attention','error','replace');
            set(AVG.hTrigger.hNumber,'string','0');
        else
            newTime=AVG.Params.currentTrigger(tmpTriggerNumber);
            AVG.Params.startTime=newTime+AVG.Params.triggerOffset;
            set(AVG.hNav.startTimeEdit,'string',num2str(AVG.Params.startTime));
            CallbackNavTimeEdit()
        end
        
    end

    function CallbackTriggerDirectionPush(hObj,Event,direction)        
        if direction==1 %forward direction
            if AVG.Params.startTime(1)<AVG.Params.currentTriggerSorted(end)+AVG.Params.triggerOffset
                tmpNumber=AVG.Params.currentTriggerSortedOrder(find(AVG.Params.currentTriggerSorted>(0.001+AVG.Params.startTime(1)-AVG.Params.triggerOffset),1,'first'));
                if isempty(tmpNumber) %this is only reached in extreme cases due to rounding problems
                    msgbox('Ending of trigger list reached','Attention','error','replace');
                    return;
                end
            else
                msgbox('Ending of trigger list reached','Attention','error','replace');
                return;
            end
        else %backward direction
            if AVG.Params.startTime(1)>(AVG.Params.currentTriggerSorted(1)+AVG.Params.triggerOffset)
                tmpNumber=AVG.Params.currentTriggerSortedOrder(find(AVG.Params.currentTriggerSorted<(AVG.Params.startTime(1)-AVG.Params.triggerOffset-0.001),1,'last'));
                if isempty(tmpNumber) %this is only reached in extreme cases due to rounding problems
                    msgbox('Beginning of trigger list reached','Attention','error','replace');
                    return;
                end
            else
                msgbox('Beginning of trigger list reached','Attention','error','replace');
                return;
            end
        end
        newTime=AVG.Params.currentTrigger(tmpNumber);
        set(AVG.hTrigger.hNumber,'string',num2str(tmpNumber));
        AVG.Params.startTime=newTime+AVG.Params.triggerOffset;
        set(AVG.hNav.startTimeEdit,'string',num2str(AVG.Params.startTime));
        CallbackNavTimeEdit()
    end

    function CallbackExportSelectedTrigger(hObj,Event)
        assignin('base','exportedTrigger',AVG.Params.currentTrigger);
    end

    function CallbackSendTriggerToStartTime(hObj,Event)
        if AVG.Params.nTriggers>0
            set(AVG.hNav.startTimeEdit,'string',num2str(AVG.Params.currentTrigger+AVG.Params.triggerOffset))
            CallbackNavTimeEdit()
        else
            msgbox('There are no selected triggers','Attention','error','replace');
        end
    end

    function CallbackDeleteTrigger(hObj,Event)
        %delete selected trigger
        if AVG.Params.nTriggers~=1
            selectedTrigger=find(cell2mat(get(AVG.hTrigger.hCheck,'value')));
        else
            selectedTrigger=1;
        end
        AVG.Params.triggers(selectedTrigger)=[];
        AVG.Params.nTriggers=AVG.Params.nTriggers-1;
        createTriggerGUI();
    end

    function CallbackAddTrigSubPopPush(hObj,Event)
        tmp=str2num(get(AVG.hTrigger.selectSubPopEdit,'string'));
        if isnumeric(tmp) && min(tmp)>0 && max(tmp)<=numel(AVG.Params.currentTrigger)
            AVG.Params.triggers=[AVG.Params.triggers {AVG.Params.currentTrigger(tmp)}];
            AVG.Params.nTriggers=AVG.Params.nTriggers+1;
            createTriggerGUI();
            CallbackTriggerCheckbox([],[],AVG.Params.nTriggers)
        else
            msgbox('Entered values are not valid triger places','Attention','error','replace');
            set(AVG.hTrigger.selectSubPopEdit,'string','')
        end
    end
    
    function CallbackTriggerCheckbox(hObj,Event,triggerNumber)
        set(AVG.hTrigger.hCheck(:),'value',0);
        set(AVG.hTrigger.hCheck(triggerNumber),'value',1);
        AVG.Params.currentTrigger=AVG.Params.triggers{triggerNumber};
        [AVG.Params.currentTriggerSorted,AVG.Params.currentTriggerSortedOrder]=sort(AVG.Params.currentTrigger);
        AVG.Params.nCurrentTriggers=numel(AVG.Params.currentTrigger);
    end

    function CallbackLoadTriggerData(hObj,Event)
        if hObj.Value
            if ~isempty(AVG.recordingObj)
                initializeTriggers;
            else
                hObj.Value=false;
                msgbox('Can not load triggers since no recording was selected! Select and try again','Attention','error','replace');
            end
        end
    end

    function CallbackPlotTriggerData(hObj,Event)
        if hObj.Value==0
            AVG.plotData.trigMarks=[];
        else
            if AVG.Params.nCurrentTriggers>0
                AVG.plotData.trigMarks=AVG.Params.currentTriggerSorted(AVG.Params.currentTriggerSorted>AVG.Params.startTime(1) & AVG.Params.currentTriggerSorted<=(AVG.Params.startTime(1)+AVG.Params.window))-AVG.Params.startTime(1);
            end
        end
        AVG.plotData.replot;
    end
%% %%%%%%%%%%%%%%%%%%%%%%%%%% Create GUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function createAVGUI(figH)
        % Construct main GUI screen
        %% Open a window and add some menus
        if nargin<1
            AVG.hMainFigure.scrsz = get(0,'ScreenSize');
            AVG.hMainFigure.hFigure = figure('Position',[AVG.hMainFigure.scrsz(3)*0.01 AVG.hMainFigure.scrsz(4)*0.07 AVG.hMainFigure.scrsz(3)*0.98 AVG.hMainFigure.scrsz(4)*0.85], ...
                'Name','Activity viewer', 'NumberTitle','off', 'MenuBar','none', 'Toolbar','none', 'HandleVisibility','off','CloseRequestFcn',@closeMainGUIFigure);
        end
        
        % define zoom options
        AVG.hMainFigure.hZoom = zoom(AVG.hMainFigure.hFigure);
        set(AVG.hMainFigure.hZoom,'Enable','on','Motion','Both','RightClickAction','PostContextMenu');
        
        % Set default panel color
        %uix.set(AVG.hMainFigure.hFigure, 'DefaultBoxPanelTitleColor',[0.7 1.0 0.7]);
        
        % set file menus
        AVG.hMenubar.hFile = uimenu(AVG.hMainFigure.hFigure, 'Label', 'File' );
        AVG.hMenubar.exportData=uimenu(AVG.hMenubar.hFile, 'Label', 'Export data', 'Callback', @CallbackExportData);
        %AVG.hMenubar.fileLoad=uimenu( AVG.hMenubar.hFile, 'Label', 'Load', 'Callback', @CallbackMenuLoad);
        %AVG.hMenubar.fileExit=uimenu( AVG.hMenubar.hFile, 'Label', 'Exit', 'Callback', @CallbackMenuExit);
        AVG.hMenubar.addSyncAV=uimenu(AVG.hMenubar.hFile, 'Label', 'Add synced AV', 'Callback', @addSyncAV);
        AVG.hMenubar.addSyncVideo=uimenu(AVG.hMenubar.hFile, 'Label', 'Add synced video', 'Callback', @addSyncVideo);

        
        % set plot meuns
        AVG.hMenubar.hPlotType = uimenu(AVG.hMainFigure.hFigure, 'Label', 'Plot' );
        for i=1:length(AVG.Params.PlotList)
            AVG.hMenubar.hPlot.(AVG.Params.PlotList{i})=uimenu(AVG.hMenubar.hPlotType,...
                'Label', AVG.Params.PlotList{i}, 'Checked','off', 'Callback', @CallbackPlotTypeSelection);
        end
        set(AVG.hMenubar.hPlot.(AVG.Params.currentPlotName), 'Checked','on');
        
        % set data type meuns
        AVG.hMenubar.hDataType = uimenu(AVG.hMainFigure.hFigure, 'Label', 'Load data' );
        for i=1:length(AVG.Params.dataRecordingClasses)
            AVG.hMenubar.hRecording.(AVG.Params.dataRecordingClasses{i})=uimenu(AVG.hMenubar.hDataType,...
                'Label', AVG.Params.dataRecordingClasses{i}, 'Checked','off', 'Callback', @CallbackDataTypeSelection);
        end
        
        AVG.hMenubar.hCredits = uimenu(AVG.hMainFigure.hFigure, 'Label', 'Other' );
        AVG.hMenubar.hMainDeveloper=uimenu(AVG.hMenubar.hCredits, 'Label', 'Developed by', 'Callback', @CallbackCredits);
        
        % Arrange the main interface windows
        AVG.hMainFigure.hMainWindow = uix.HBox('Parent',AVG.hMainFigure.hFigure, 'Spacing',4);
        AVG.hMainFigure.hLeftBox = uix.VBoxFlex('Parent',AVG.hMainFigure.hMainWindow, 'Spacing',4, 'Padding',4);
        AVG.hMainFigure.hMidGrid = uix.Grid('Parent',AVG.hMainFigure.hMainWindow, 'Spacing',4, 'Padding',7);
        AVG.hMainFigure.hRightBox = uix.VBoxFlex('Parent',AVG.hMainFigure.hMainWindow, 'Spacing',4, 'Padding',4);
        set(AVG.hMainFigure.hMainWindow, 'Widths',[-2 -8 -2]);
        
        % Set left box
        AVG.hGen.General = uix.Panel('Parent',AVG.hMainFigure.hLeftBox, 'Title','General');
        AVG.hNav.MainPanel = uix.Panel('Parent',AVG.hMainFigure.hLeftBox, 'Title','Navigation');
        AVG.hCh.ChannelSelectionPanel = uix.Panel('Parent',AVG.hMainFigure.hLeftBox, 'Title','Channel selection');
        AVG.hPlotProp.hMainPanel = uix.Panel('Parent',AVG.hMainFigure.hLeftBox, 'Title','Plot specific properties');
        set(AVG.hMainFigure.hLeftBox, 'Heights',[-1 -2 -2 -3]);
        
        % Set mid box - main plot
        AVG.hMainFigure.dummyEmptyBoxes(1)=uix.Empty('Parent', AVG.hMainFigure.hMidGrid,'Visible','off');
        AVG.hMainFigure.dummyEmptyBoxes(2)=uix.Empty('Parent', AVG.hMainFigure.hMidGrid,'Visible','off');
        AVG.hMainFigure.dummyEmptyBoxes(3)=uix.Empty('Parent', AVG.hMainFigure.hMidGrid,'Visible','off');
        AVG.hMainFigure.hTimeSlider=uicontrol('Parent',AVG.hMainFigure.hMidGrid, 'Style','slider', 'Min',0, 'Max',1, 'Value',0, 'Callback',@CallbackTimeSlider);
        AVG.hMainFigure.hMainAxis=axes('Parent', AVG.hMainFigure.hMidGrid, 'ActivePositionProperty','Position', 'FontSize',8,'NextPlot','add','ColorOrderIndex',1);
        AVG.hMainFigure.dummyEmptyBoxes(4)=uix.Empty('Parent', AVG.hMainFigure.hMidGrid,'Visible','off');
        set(AVG.hMainFigure.hMidGrid,'Widths',[45 -1],'Heights', [AVG.Params.scrollBarPixels -1 35]);
        
        % Set right box
        AVG.hProcessing.hMainPanel = uix.TabPanel('Parent',AVG.hMainFigure.hRightBox, 'Padding',2, 'TabWidth',90);
        AVG.hFilter.hMainPanel=uix.Panel('Parent',AVG.hProcessing.hMainPanel, 'Padding',4);
        AVG.hAnalysis.hMainPanel=uix.Panel('Parent',AVG.hProcessing.hMainPanel, 'Padding',4);
        set(AVG.hProcessing.hMainPanel,'TabTitles',{'Filter','Analysis'}, 'Selection',2);
        
        AVG.hTrigger.hMainTriggerPanel = uix.Panel('Parent',AVG.hMainFigure.hRightBox, 'Title','Triggers');
        set(AVG.hMainFigure.hRightBox, 'Heights',[-2,-3]);
        %% construch general GUI
        AVG.hGen.mainVBox=uix.VBox('Parent', AVG.hGen.General, 'Padding', 2, 'Spacing', 5);
        AVG.hGen.messageBox=uicontrol('Parent', AVG.hGen.mainVBox, 'HorizontalAlignment','left', 'Style','text',...
            'String','Ready','ForegroundColor','k','FontSize',12,'FontWeight','Bold');
        
        AVG.hGen.generalGrid=uix.Grid('Parent', AVG.hGen.mainVBox, 'Padding', 2, 'Spacing', 5);
        AVG.hGen.refreshPlotPush=uicontrol('Parent', AVG.hGen.generalGrid, 'Callback',{@CallbackRefreshPlotPush}, 'Style','push', 'String','Refresh plot');
        AVG.hGen.restartPlotPush=uicontrol('Parent', AVG.hGen.generalGrid, 'Callback',{@CallbackRestartPlotPush}, 'Style','push', 'String','Restart plot');
        AVG.hGen.exportPlotPush=uicontrol('Parent', AVG.hGen.generalGrid, 'Callback',{@CallbackExportPlotPush}, 'Style','push', 'String','Export plot');
                
        set(AVG.hGen.generalGrid, 'Widths',-1,'Heights', -1 );
        set(AVG.hGen.mainVBox,'Heights',[-1,-3]);
        %% Construct Navigation GUI
        AVG.hNav.mainBox=uix.VBox('Parent', AVG.hNav.MainPanel, 'Padding', 5, 'Spacing', 10);
        
        AVG.hNav.timeControlBox=uix.HBox('Parent', AVG.hNav.mainBox, 'Padding', 5, 'Spacing', 10);
        AVG.hNav.hFastBackPush=uicontrol('Parent', AVG.hNav.timeControlBox, 'Callback',{@CallbackNavTimeFrameMovePush,-1}, 'Style','push', 'String','<<','FontSize',12,'FontWeight','Bold','BackgroundColor',[0.6 0.6 1]);
        AVG.hNav.hBackPush=uicontrol('Parent', AVG.hNav.timeControlBox, 'Callback',{@CallbackNavTimeFrameMovePush,-AVG.Params.shortTimeFrameMovement}, 'Style','push', 'String','<','FontSize',12,'FontWeight','Bold','BackgroundColor',[1 0.6 0.6]);
        AVG.hNav.ZoomBox=uix.VBox('Parent', AVG.hNav.timeControlBox, 'Padding', 5, 'Spacing', 10);
        AVG.hNav.hZoomInPush=uicontrol('Parent', AVG.hNav.ZoomBox, 'Callback',{@CallbackNavZoomPush,0.5}, 'Style','push','FontSize',12,'FontWeight','Bold','String','+');
        AVG.hNav.hZoomOutPush=uicontrol('Parent', AVG.hNav.ZoomBox, 'Callback',{@CallbackNavZoomPush,2}, 'Style','push','FontSize',12,'FontWeight','Bold','String','-');
        AVG.hNav.hForwardPush=uicontrol('Parent', AVG.hNav.timeControlBox, 'Callback',{@CallbackNavTimeFrameMovePush,AVG.Params.shortTimeFrameMovement}, 'Style','push', 'String','>','FontSize',12,'FontWeight','Bold','BackgroundColor',[1 0.6 0.6]);
        AVG.hNav.hFastForwardPush=uicontrol('Parent', AVG.hNav.timeControlBox, 'Callback',{@CallbackNavTimeFrameMovePush,1}, 'Style','push', 'String','>>','FontSize',12,'FontWeight','Bold','BackgroundColor',[0.6 0.6 1]);
        
        set(AVG.hNav.timeControlBox, 'Widths',[-3,-3,-4,-3,-3]);
        
        AVG.hNav.timeControlGrid=uix.Grid('Parent', AVG.hNav.mainBox, 'Padding', 5, 'Spacing', 10);
        
        AVG.hNav.TimeTxt=uicontrol('Parent', AVG.hNav.timeControlGrid, 'HorizontalAlignment','left', 'Style','text', 'String','start time [ms]');
        AVG.hNav.WindowTxt=uicontrol('Parent', AVG.hNav.timeControlGrid, 'HorizontalAlignment','left','Style','text', 'String','window [ms]');
        
        AVG.hNav.startTimeEdit=uicontrol('Parent', AVG.hNav.timeControlGrid, 'Callback',{@CallbackNavTimeEdit},'Style','edit', 'String',AVG.Params.defaultStartTime);
        AVG.hNav.WindowEdit=uicontrol('Parent', AVG.hNav.timeControlGrid, 'Callback',{@CallbackNavWindowEdit}, 'Style','edit', 'String',AVG.Params.defaultWindow);

        set(AVG.hNav.timeControlGrid, 'Widths',[-2 -1],'Heights', [-1 -1] );
        
        set(AVG.hNav.mainBox, 'Heights',[-3,-2]);
        %% Construct channel GUI
        AVG.hCh.mainBox=uix.VBox('Parent', AVG.hCh.ChannelSelectionPanel, 'Padding', 5, 'Spacing', 10);
        
        AVG.hCh.AddChannelsEdit=uicontrol('Parent', AVG.hCh.mainBox, 'Callback',{@CallbackChAddChannelsEdit},'Style','edit', 'String','Add channels');
        AVG.hCh.removeChannelsEdit=uicontrol('Parent', AVG.hCh.mainBox, 'Callback',{@CallbackChRemoveChannelsEdit}, 'Style','edit', 'String','Remove channels');
        AVG.hCh.keepOnlyChannelsEdit=uicontrol('Parent', AVG.hCh.mainBox, 'Callback',{@CallbackChKeepOnlyChannelsEdit}, 'Style','edit', 'String','Keep only channels');
        AVG.hCh.enableAllChannelsPush=uicontrol('Parent', AVG.hCh.mainBox, 'Callback',{@CallbackChEnableAllChannelsPush}, 'Style','push', 'String','Enable all channels');
        AVG.hCh.manualChannelSelectPush=uicontrol('Parent', AVG.hCh.mainBox, 'Callback',{@CallbackChManualChannelSelectPush}, 'Style','push', 'String','Select manually');
        
        %AVG.hCh.empty1=uix.Empty('Parent', AVG.hCh.mainBox);
        
        set(AVG.hCh.mainBox, 'Heights',[-1,-1,-1,-1,-1]);
        %% Construct filter GUI
        % set filter panel
        %fileds are arranged as a N x 2 grid of text and input
        
        AVG.hFilter.hMainVBox= uix.VBox( 'Parent', AVG.hFilter.hMainPanel, 'Spacing',4, 'Padding',4);
        AVG.hFilter.hMainGrid= uix.Grid( 'Parent', AVG.hFilter.hMainVBox, 'Spacing',10 );
        
        %text column
        AVG.hFilter.hFilterTxt=uicontrol('Parent', AVG.hFilter.hMainGrid, 'HorizontalAlignment','left', 'Style','text', 'String','Filter');
        AVG.hFilter.hDesignTxt=uicontrol('Parent', AVG.hFilter.hMainGrid, 'HorizontalAlignment','left', 'Style','text', 'String','Filter design');
        AVG.hFilter.hOrderTxt=uicontrol('Parent', AVG.hFilter.hMainGrid, 'HorizontalAlignment','left','Style','text', 'String','Filter order');
        AVG.hFilter.hSamplingFrequencyTxt=uicontrol('Parent', AVG.hFilter.hMainGrid, 'HorizontalAlignment','left','Style','text', 'String','Sampling Freq. [Hz]');
        AVG.hFilter.hLowPassStopTxt=uicontrol('Parent', AVG.hFilter.hMainGrid, 'HorizontalAlignment','left','Style','text', 'String','Lowpass stop Freq. [Hz]');
        AVG.hFilter.hLowPassPassTxt=uicontrol('Parent', AVG.hFilter.hMainGrid, 'HorizontalAlignment','left','Style','text', 'String','Lowpass pass Freq. [Hz]');
        AVG.hFilter.hAttenuationLowPassTxt=uicontrol('Parent', AVG.hFilter.hMainGrid, 'HorizontalAlignment','left','Style','text', 'String','Lowpass attenuation [dB]');
        AVG.hFilter.hHighPassPassTxt=uicontrol('Parent', AVG.hFilter.hMainGrid, 'HorizontalAlignment','left','Style','text', 'String','Highpass pass Freq. [Hz]');
        AVG.hFilter.hHighPassStopTxt=uicontrol('Parent', AVG.hFilter.hMainGrid, 'HorizontalAlignment','left','Style','text', 'String','Highpass stop Freq. [Hz]');
        AVG.hFilter.hAttenuationHighTxt=uicontrol('Parent', AVG.hFilter.hMainGrid, 'HorizontalAlignment','left','Style','text', 'String','High attenuation [dB]');
        AVG.hFilter.hRippleTxt=uicontrol('Parent', AVG.hFilter.hMainGrid, 'HorizontalAlignment','left','Style','text', 'String','Ripple (band) [dB]');
        AVG.hFilter.hSamplingRatio=uicontrol('Parent', AVG.hFilter.hMainGrid, 'HorizontalAlignment','left','Style','text', 'String','Sampling ratio');
        AVG.hFilter.hPaddingSamplesTxt=uicontrol('Parent', AVG.hFilter.hMainGrid, 'HorizontalAlignment','left','Style','text', 'String','Padding [ms]');
        
        %input column
        newPopUpMenu=[AVG.Params.filterNames];
        AVG.hFilter.hMethodPopup=uicontrol('Parent', AVG.hFilter.hMainGrid, 'Callback',{@CallbackFilterMethodPopup},'Style','popup', 'String',newPopUpMenu,'Min',1,'Max',numel(newPopUpMenu));
        AVG.hFilter.hDesignPopup=uicontrol('Parent', AVG.hFilter.hMainGrid, 'Callback',{@CallbackFilterDesignPopup},'Style','popup', 'String',AVG.Params.filterDesigns,'value',AVG.Params.defaultFilterValue);
        AVG.hFilter.hOrderEdit=uicontrol('Parent', AVG.hFilter.hMainGrid, 'Callback',{@CallbackFilterOrderEdit}, 'Style','edit', 'String',AVG.Params.defaultFilterOrder);
        AVG.hFilter.hSamplingFrequencyEdit=uicontrol('Parent', AVG.hFilter.hMainGrid, 'Callback',{@CallbackFilterSamplingFrequencyEdit}, 'Style','edit', 'String',num2str(AVG.Params.samplingFrequency(1)));
        AVG.hFilter.hLowPassStopEdit=uicontrol('Parent', AVG.hFilter.hMainGrid, 'Callback',{@CallbackFilterLowPassStopEdit}, 'Style','edit', 'String',AVG.Params.defaultFilterLowPassStopFreq);
        AVG.hFilter.hLowPassPassEdit=uicontrol('Parent', AVG.hFilter.hMainGrid, 'Callback',{@CallbackFilterLowPassPassEdit}, 'Style','edit', 'String',AVG.Params.defaultFilterLowPassPassFreq);
        AVG.hFilter.hLowpassAttenuationEdit=uicontrol('Parent', AVG.hFilter.hMainGrid, 'Callback',{@CallbackFilterLowPassAttenuationEdit}, 'Style','edit', 'String',AVG.Params.defaultAttenuationInLowpass);
        AVG.hFilter.hHighPassPassEdit=uicontrol('Parent', AVG.hFilter.hMainGrid, 'Callback',{@CallbackFilterHighPassPassEdit}, 'Style','edit', 'String',AVG.Params.defaultFilterHighPassPassFreq);
        AVG.hFilter.hHighPassStopEdit=uicontrol('Parent', AVG.hFilter.hMainGrid, 'Callback',{@CallbackFilterHighPassStopEdit}, 'Style','edit', 'String',AVG.Params.defaultFilterHighPassStopFreq);
        AVG.hFilter.hHighpassAttenuationEdit=uicontrol('Parent', AVG.hFilter.hMainGrid, 'Callback',{@CallbackFilterHighPassAttenuationEdit}, 'Style','edit', 'String',AVG.Params.defaultAttenuationInHighpass);
        AVG.hFilter.hRippleEdit=uicontrol('Parent', AVG.hFilter.hMainGrid, 'Callback',{@CallbackFilterRippleEdit}, 'Style','edit', 'String',AVG.Params.defaultRippleInPassband);
        
        AVG.hFilter.hSamplingRatioEdit=uicontrol('Parent', AVG.hFilter.hMainGrid, 'Callback',{@CallbackFilterSamplingRatioEdit}, 'Style','edit', 'String',AVG.Params.defaultFilterSamplingRatio);
        AVG.hFilter.hPaddingSamplesEdit=uicontrol('Parent', AVG.hFilter.hMainGrid, 'Callback',{@CallbackPaddingSamplesEdit}, 'Style','edit', 'String',AVG.Params.defaultPadding);
        
        %arrange columns into a grid
        set(AVG.hFilter.hMainGrid, 'Widths',[-2 -1],'Heights', 19*ones(1,13) );
        
        AVG.hFilter.hPlotFilterPush=uicontrol('Parent', AVG.hFilter.hMainVBox, 'Callback',{@CallbackPlotFilterPush}, 'Style','push', 'String','Plot filter');
        set(AVG.hFilter.hMainVBox, 'Heights',[-1 30]);
        %% Construct analysis GUI
        AVG.hAnalysis.hMainVbox= uix.VBox( 'Parent', AVG.hAnalysis.hMainPanel, 'Spacing',5 );
        AVG.hAnalysis.hMainGrid= uix.Grid( 'Parent', AVG.hAnalysis.hMainVbox, 'Spacing',10 );        
        for i=1:numel(AVG.Params.analysisMethodsNames)
            AVG.hAnalysis.hToggleButtons(i)=uicontrol('Parent',AVG.hAnalysis.hMainGrid,'Callback',{@CallbackAnaysisTogglePush,i},...
               'Style','toggle','String',[AVG.Params.analysisMethodsNames{i} ' (-)']);
        end
        set(AVG.hAnalysis.hMainGrid,'Widths',[-1 -1]);
        AVG.hAnalysis.hRemoveAll=uicontrol('Parent',AVG.hAnalysis.hMainVbox,'Callback',{@CallbackRemoveAllAnalysisPush}, 'Style','push', 'String','Remove all','BackgroundColor',[0.5 0.7 0.5]);
        set(AVG.hAnalysis.hMainVbox,'Heights',[-1 40]);
        %% Construct Trigger Box
        AVG.hTrigger.MainVBox=uix.VBox('Parent', AVG.hTrigger.hMainTriggerPanel, 'Padding', 4, 'Spacing', 4);
        
        AVG.hTrigger.navigationHBox=uix.HBox('Parent', AVG.hTrigger.MainVBox, 'Padding', 4, 'Spacing', 4);
        AVG.hTrigger.hForward=uicontrol('Parent',AVG.hTrigger.navigationHBox,'Callback',{@CallbackTriggerDirectionPush,-1}, 'Style','push', 'String','<<');
        AVG.hTrigger.hNumber=uicontrol('Parent',AVG.hTrigger.navigationHBox,'Callback',@CallbackTriggerNumberEdit, 'Style','edit', 'String','0');
        AVG.hTrigger.hBackward=uicontrol('Parent',AVG.hTrigger.navigationHBox,'Callback',{@CallbackTriggerDirectionPush,1}, 'Style','push', 'String','>>');
        set(AVG.hTrigger.navigationHBox, 'Widths',[-1 40 -1]);
        
        AVG.hTrigger.hLoadAndPlotTrigHBox=uix.HBox('Parent', AVG.hTrigger.MainVBox, 'Padding', 4, 'Spacing', 4);
        AVG.hTrigger.hGetTrigFromRawFiles=uicontrol('Parent', AVG.hTrigger.hLoadAndPlotTrigHBox,'Callback',{@CallbackLoadTriggerData}, 'HorizontalAlignment','left','Style','push', 'String','load trig.','value',AVG.Params.loadTriggerDefault,'Tooltip','load triggers from recording file');
        AVG.hTrigger.hPlotTrigs=uicontrol('Parent', AVG.hTrigger.hLoadAndPlotTrigHBox,'Callback',@CallbackPlotTriggerData, 'HorizontalAlignment','left','Style','check', 'String','plot trig.','value',0,'Tooltip','plot triggers on traces (works only on activity trace plot)');
        
        AVG.hTrigger.OffsetHBox=uix.HBox('Parent', AVG.hTrigger.MainVBox, 'Padding', 4, 'Spacing', 4);
        AVG.hTrigger.hOffsetTxt=uicontrol('Parent', AVG.hTrigger.OffsetHBox, 'HorizontalAlignment','left','Style','text', 'String','Offset [ms]');
        AVG.hTrigger.hOffsetEdit=uicontrol('Parent', AVG.hTrigger.OffsetHBox, 'Style','edit', 'String',0,'Callback',@CallbackTrigOffset);
        set(AVG.hTrigger.OffsetHBox, 'Widths',[-1 -2]);
        
        AVG.hTrigger.selectSubPopulationHBox=uix.HBox('Parent', AVG.hTrigger.MainVBox, 'Padding', 4, 'Spacing', 4);
        AVG.hTrigger.selectSubPopPush=uicontrol('Parent', AVG.hTrigger.selectSubPopulationHBox,'String','Add sub. pop.','Callback',@CallbackAddTrigSubPopPush,'HorizontalAlignment','left','Style','push');
        AVG.hTrigger.selectSubPopEdit=uicontrol('Parent', AVG.hTrigger.selectSubPopulationHBox,'Style','edit', 'String','');
        set(AVG.hTrigger.selectSubPopulationHBox, 'Widths',[-1 -1]);
        
        AVG.hTrigger.manualSetHBox=uix.HBox('Parent', AVG.hTrigger.MainVBox, 'Padding', 4, 'Spacing', 4);
        AVG.hTrigger.manualSetPush=uicontrol('Parent', AVG.hTrigger.manualSetHBox, 'Style','push', 'String','Add manually','Callback',@CallbackManualSetPush);
        AVG.hTrigger.manualSetEdit=uicontrol('Parent', AVG.hTrigger.manualSetHBox, 'Style','edit', 'String','');
        set(AVG.hTrigger.manualSetHBox, 'Widths',[-1 -2]);
        
        AVG.hTrigger.manualLoadHBox=uix.HBox('Parent', AVG.hTrigger.MainVBox, 'Padding', 4, 'Spacing', 4);
        AVG.hTrigger.manualLoadTxtTxt=uicontrol('Parent', AVG.hTrigger.manualLoadHBox, 'HorizontalAlignment','left','Style','text', 'String','Load trigger: ');
        AVG.hTrigger.manualLoadFile=uicontrol('Parent', AVG.hTrigger.manualLoadHBox,'Style','push','String','file','Callback',@CallbackManualLoadFile,'Tooltip','load triggers from matlab file with a saved cell array');
        AVG.hTrigger.manualLoadVariable=uicontrol('Parent', AVG.hTrigger.manualLoadHBox,'Style','push','String','var','Callback',@CallbackManualLoadVariable,'Tooltip','load triggers from cell array variable on workspace');
        AVG.hTrigger.manualLoadSpikes=uicontrol('Parent', AVG.hTrigger.manualLoadHBox,'Style','push','String','t-ic','Callback',@CallbackManualLoadSpikes,'Tooltip','load triggers spike sorting mat file with t-ic format');

        AVG.hTrigger.exportTriggerHBox=uix.HBox('Parent', AVG.hTrigger.MainVBox, 'Padding', 4, 'Spacing', 4);
        AVG.hTrigger.exportSelectedTrigger=uicontrol('Parent', AVG.hTrigger.exportTriggerHBox, 'Style','push', 'String','Export selected trigger','Callback',@CallbackExportSelectedTrigger);
        AVG.hTrigger.sendTriggerToStartTimePush=uicontrol('Parent', AVG.hTrigger.exportTriggerHBox, 'Style','push', 'String','send to start times','Callback',@CallbackSendTriggerToStartTime,'ForegroundColor','r');
        
        AVG.hTrigger.manualLoadFile=uicontrol('Parent',AVG.hTrigger.MainVBox,'Style','push','String','delete trigger','Callback',@CallbackDeleteTrigger);
        
        %delete this handle for terminating the trigger GUI and replacing it with a new one
        AVG.hTrigger.MainGrid=uix.Grid('Parent',AVG.hTrigger.MainVBox, 'Spacing',4, 'Padding',4);
        
        set(AVG.hTrigger.MainVBox, 'Heights',[50 25 25 25 25 25 50 25 -1]);
        %this line should be after the startTimeEdit uicontrol is defined! 
        
        %mouse_figure(AVG.hMainFigure.hFigure);
    end %createGUI
%% %%%%%%%%%%%%%%%%%%%%%%%%%% Create Trigger GUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function createTriggerGUI() %the function sets the first trigger to be default
        if isfield(AVG.hTrigger,'MainGrid') %delete trigger graphical objects
            delete(AVG.hTrigger.MainGrid);
            AVG.hTrigger.MainGrid=uix.Grid('Parent',AVG.hTrigger.MainVBox, 'Spacing',4, 'Padding',4);
        end
        AVG.hTrigger.hCheck=[];
        %plot graphical checkboxes for active trigger
        for i=1:AVG.Params.nTriggers
            AVG.hTrigger.hCheck(i)=uicontrol('Parent',AVG.hTrigger.MainGrid,'Style','checkbox',...
                'String',[num2str(i) ' - N=' num2str(numel(AVG.Params.triggers{i}))],...
                'value',i==1,'Callback',{@CallbackTriggerCheckbox,i});
        end
        if AVG.Params.nTriggers>0
            set(AVG.hTrigger.MainGrid, 'Widths', -1*ones(1,min(3,ceil(AVG.Params.nTriggers/AVG.Params.maxNumberOfTriggerInColumn))), 'Heights', 20*ones(1,min(AVG.Params.maxNumberOfTriggerInColumn,AVG.Params.nTriggers)));
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create Sync video GUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function createSyncVideoGUI()
        createIcons;
        AVG.hVideoSyncFigure.hFigure = figure('Position',[AVG.hMainFigure.scrsz(3)*0.5 AVG.hMainFigure.scrsz(4)*0.4 AVG.hMainFigure.scrsz(3)*0.48 AVG.hMainFigure.scrsz(4)*0.45], ...
            'Name','Activity viewer - video sync', 'NumberTitle','off', 'MenuBar','none', 'Toolbar','figure', 'HandleVisibility','off','CloseRequestFcn',@closeSyncedVideoFigure);
        %To make the zoom toolbar not visible change 'Toolbar' to 'none'

        % Arrange the main interface windows
        AVG.hVideoSyncFigure.hMainVBox = uix.VBox('Parent',AVG.hVideoSyncFigure.hFigure, 'Spacing',3, 'Padding',1);

        AVG.hVideoSyncFigure.hButtonPanelTop = uix.Panel('Parent',AVG.hVideoSyncFigure.hMainVBox, 'Title','Controls');
        AVG.hVideoSyncFigure.hVideoPanel = uix.Panel('Parent',AVG.hVideoSyncFigure.hMainVBox, 'Title','Video');
        set(AVG.hVideoSyncFigure.hMainVBox, 'Heights',[-2 -10]);

        AVG.hVideoSyncFigure.hVideoAxis=axes('Parent', AVG.hVideoSyncFigure.hVideoPanel,'Position',[0.025 0.025 0.95 0.95]);
        %AVG.hVideoSyncFigure.hVideoAxis.Toolbar.Visible = 'on';
        %hold(AVG.hVideoSyncFigure.hVideoAxis,'off');

        %High horizontalPanel
        AVG.hVideoSyncFigure.hControlPanelVBox = uix.VBox('Parent',AVG.hVideoSyncFigure.hButtonPanelTop);
        AVG.hVideoSyncFigure.hTopButtonHBox = uix.HBox('Parent',AVG.hVideoSyncFigure.hControlPanelVBox, 'Spacing',4, 'Padding',2);
        AVG.hVideoSyncFigure.hBottomButtonHBox = uix.HBox('Parent',AVG.hVideoSyncFigure.hControlPanelVBox, 'Spacing',4, 'Padding',2);
        set(AVG.hVideoSyncFigure.hControlPanelVBox, 'Heights',[-1 -1]);

        AVG.hVideoSyncFigure.hLoadVideoPush=uicontrol('Parent', AVG.hVideoSyncFigure.hTopButtonHBox, 'Callback',@CallbackLoadVideoPush, 'Style','push',...
            'Tooltip','Load video','CData',AVG.Params.Icons.load,'FontSize',12,'FontWeight','Bold','BackgroundColor',[0.8 0.8 0.8]);
        AVG.hVideoSyncFigure.hLoadVideoEdit=uicontrol('Parent', AVG.hVideoSyncFigure.hTopButtonHBox,'Style','edit',...
            'String','videoFileName','Callback',@CallbackLoadVideoEdit);
        AVG.hVideoSyncFigure.hCheckSyncPush=uicontrol('Parent', AVG.hVideoSyncFigure.hTopButtonHBox, 'Callback',@CallbackCheckSyncPush, 'Style','push',...
            'Tooltip','Check sync.','CData',AVG.Params.Icons.syncIcon,'FontSize',12,'FontWeight','Bold','BackgroundColor',[0.8 0.8 0.8]);
        AVG.hVideoSyncFigure.hValidSyncTriggersEdit=uicontrol('Parent', AVG.hVideoSyncFigure.hTopButtonHBox, 'Callback',@CallbackValidSyncTriggersEdit, ...
            'Style','edit', 'String','0-for-uniform','FontSize',12,'FontWeight','Bold','BackgroundColor',[0.8 0.8 0.8],'Tooltip','If 0 is entered, I will assume that frames are lost at a uniform rate');
        AVG.hVideoSyncFigure.hFrameRateEdit=uicontrol('Parent', AVG.hVideoSyncFigure.hTopButtonHBox, ...
            'Callback',@CallbackFrameRateEdit, 'Style','edit', 'String','Frame rate','FontSize',12,'FontWeight','Bold','BackgroundColor',[0.8 0.8 0.8]);
        AVG.hVideoSyncFigure.hPlaySpeedupPush=uicontrol('Parent', AVG.hVideoSyncFigure.hTopButtonHBox, 'Callback',@CallbackPlaySpeedupPush, 'Style','push',...
            'String','X1','Tooltip','Play video faster (push to switch)','FontSize',10,'FontWeight','Bold','BackgroundColor',[0.8 0.8 0.8]);
        AVG.hVideoSyncFigure.hPlayVideoPush=uicontrol('Parent', AVG.hVideoSyncFigure.hTopButtonHBox, 'Callback',@CallbackRunVideoPush, 'Style','push',...
            'Tooltip','Run video','CData',AVG.Params.Icons.playForwardIcon,'FontSize',12,'FontWeight','Bold','BackgroundColor',[0.8 0.8 0.8]);
        AVG.hVideoSyncFigure.hExportIMPlayPush=uicontrol('Parent', AVG.hVideoSyncFigure.hTopButtonHBox, 'Callback',@CallbackExportIMPlayPush, 'Style','push',...
            'Tooltip','Export video to IMPlay','String','Exp.','FontSize',12,'FontWeight','Bold','BackgroundColor',[0.8 0.8 0.8]);
        set(AVG.hVideoSyncFigure.hTopButtonHBox, 'Widths',[-1 -6 -1 -3 -1 -1 -1 -1]);

        %Low horizontalPanel
        AVG.hVideoSyncFigure.hConvertVideo2EphysEdit=uicontrol('Parent', AVG.hVideoSyncFigure.hBottomButtonHBox, 'Callback',@CallbackConvertVideo2EphysEdit, ...
            'Style','edit', 'String','hh:mm:ss','FontSize',12,'FontWeight','Bold','BackgroundColor',[0.8 0.8 0.8],'Tooltip','Write time in video (hh:mm:ss), press enter to get time in activity viewer (ms)');
        AVG.hVideoSyncFigure.hConvertEphys2VideoEdit=uicontrol('Parent', AVG.hVideoSyncFigure.hBottomButtonHBox, 'Callback',@CallbackConvertEphys2VideoEdit, ...
            'Style','edit', 'String','time [ms]','FontSize',12,'FontWeight','Bold','BackgroundColor',[0.8 0.8 0.8],'Tooltip','Write time in ephys (ms), press enter to get time in video (hh:mm:ss)');
        set(AVG.hVideoSyncFigure.hBottomButtonHBox, 'Widths',[100 100]);
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%% Create Icons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function createIcons() %the function sets the first trigger to be default
        
        AVG.Params.Icons.playForwardIcon =...
           [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
            NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN;
            NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN;
            NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN;
            NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN;
            NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
            NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN;
            NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
            NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN;
            NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN;
            NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN;
            NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN;
            NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
        AVG.Params.Icons.playForwardIcon = repmat(AVG.Params.Icons.playForwardIcon,[1 1 3]);

        AVG.Params.Icons.playForwardFastIcon = repmat(AVG.Params.Icons.playForwardIcon,[1,2,1]);

        AVG.Params.Icons.playBackwardIcon = flip(AVG.Params.Icons.playForwardIcon,2);

        AVG.Params.Icons.playBackwardFastIcon = repmat(AVG.Params.Icons.playBackwardIcon,[1,2,1]);

        AVG.Params.Icons.pauseIcon = ...
            [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
            NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
            NaN,0  ,0  ,NaN,0  ,0  ,NaN,NaN;
            NaN,0  ,0  ,NaN,0  ,0  ,NaN,NaN;
            NaN,0  ,0  ,NaN,0  ,0  ,NaN,NaN;
            NaN,0  ,0  ,NaN,0  ,0  ,NaN,NaN;
            NaN,0  ,0  ,NaN,0  ,0  ,NaN,NaN;
            NaN,0  ,0  ,NaN,0  ,0  ,NaN,NaN;
            NaN,0  ,0  ,NaN,0  ,0  ,NaN,NaN;
            NaN,0  ,0  ,NaN,0  ,0  ,NaN,NaN;
            NaN,0  ,0  ,NaN,0  ,0  ,NaN,NaN;
            NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
            NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
        AVG.Params.Icons.pauseIcon = repmat(AVG.Params.Icons.pauseIcon,[1 1 3]);


        AVG.Params.Icons.stopIcon = ...
            [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
            NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
            NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
            NaN NaN   0   0   0   0   0   0   0 NaN NaN;
            NaN NaN   0   0   0   0   0   0   0 NaN NaN;
            NaN NaN   0   0   0   0   0   0   0 NaN NaN;
            NaN NaN   0   0   0   0   0   0   0 NaN NaN;
            NaN NaN   0   0   0   0   0   0   0 NaN NaN;
            NaN NaN   0   0   0   0   0   0   0 NaN NaN;
            NaN NaN   0   0   0   0   0   0   0 NaN NaN;
            NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
            NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
            NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
        AVG.Params.Icons.stopIcon = repmat(AVG.Params.Icons.stopIcon,[1 1 3]);

        AVG.Params.Icons.zoomIcon(:,:,1) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,254,252,252,254,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,229,167,145,129,128,150,213,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,217,151,143,153,152,137,113,112,180,249,255,255,255,255,255,255,255,255,255;255,255,255,255,225,141,176,240,250,253,236,179,116,91,198,255,255,255,255,255,255,255,255,255;255,255,255,240,152,162,253,255,249,232,220,221,182,98,111,219,255,255,255,255,255,255,255,255;255,255,255,231,127,197,255,247,231,217,206,193,172,125,87,175,255,255,255,255,255,255,255,255;255,255,255,225,109,215,254,228,216,206,194,181,172,133,87,162,255,255,255,255,255,255,255,255;255,255,255,223,99,203,242,217,207,195,189,202,206,132,87,159,255,255,255,255,255,255,255,255;255,255,255,225,104,159,212,205,190,187,213,249,204,113,87,168,255,255,255,255,255,255,255,255;253,239,239,227,129,107,177,203,160,187,224,203,136,81,108,207,255,254,255,255,255,255,255,255;255,243,239,242,208,101,103,136,136,156,160,115,72,71,85,169,252,255,255,255,255,255,255,255;255,248,239,239,250,180,107,87,83,81,88,89,79,73,83,96,184,255,255,255,255,255,255,255;255,252,239,239,248,254,193,121,94,91,106,129,147,100,87,87,96,183,255,255,255,255,255,255;255,255,242,239,244,255,255,236,231,237,245,244,207,135,89,87,87,93,189,240,239,242,249,255;255,255,246,239,240,255,252,235,227,233,237,233,223,208,139,99,87,78,86,182,242,242,250,255;255,255,250,239,239,252,248,239,239,252,252,250,235,237,225,153,98,71,71,80,197,255,254,255;255,255,254,239,239,248,244,239,246,255,255,255,239,239,247,229,147,88,71,71,88,244,254,255;255,255,255,244,239,246,242,239,249,255,255,255,239,239,248,255,231,155,88,76,108,231,245,255;255,255,255,248,239,242,239,239,254,255,255,255,239,239,248,255,255,221,156,144,188,225,236,255;255,255,255,253,239,239,239,244,255,255,255,255,239,239,248,255,252,230,195,190,216,228,244,255;255,255,255,255,243,239,239,249,255,255,255,255,239,239,248,255,254,246,242,239,243,248,254,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
        AVG.Params.Icons.zoomIcon(:,:,2) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,253,251,253,255;255,255,255,255,255,255,255,255,254,252,252,254,255,255,255,255,255,255,255,255,248,247,254,255;255,249,248,252,248,247,226,165,138,122,125,150,213,255,251,249,247,248,250,248,247,247,247,252;255,247,247,247,249,210,143,136,147,144,131,113,112,173,242,249,249,247,250,249,247,248,250,254;254,247,247,254,223,134,169,239,249,245,229,179,111,84,193,255,250,247,251,251,247,251,255,255;252,247,250,240,148,154,247,255,246,225,214,219,175,90,110,219,248,247,253,250,247,252,255,255;251,247,251,231,122,190,251,247,227,210,201,191,165,118,87,173,247,247,255,248,247,253,255,255;249,247,253,225,104,208,252,228,211,198,192,180,164,126,81,155,247,247,255,247,247,248,253,255;249,247,255,223,93,194,240,217,202,187,188,202,202,126,81,157,248,248,255,252,248,248,255,255;255,255,255,225,104,160,211,205,190,187,213,249,204,113,87,168,255,255,255,255,255,255,255,255;254,249,249,233,129,107,176,202,167,196,231,203,139,86,109,207,255,254,255,255,255,255,255,255;255,250,249,250,208,101,103,137,145,165,165,115,81,81,90,169,252,255,255,255,255,255,255,255;255,252,249,249,252,180,107,90,92,90,90,89,86,82,85,96,184,255,255,255,255,255,255,255;255,254,249,249,252,254,193,125,103,99,107,129,146,99,87,87,96,183,255,255,255,255,255,255;255,255,250,249,251,255,255,245,240,242,245,244,216,142,93,87,87,93,194,249,249,250,253,255;255,255,251,249,249,255,254,245,238,236,235,232,232,217,143,99,87,83,95,191,249,250,253,255;255,255,253,249,249,254,252,249,248,252,251,250,244,246,229,154,98,81,81,85,197,255,255,255;255,255,255,249,249,252,251,249,251,255,255,255,249,249,251,230,146,93,81,81,97,248,255,255;255,255,255,251,249,251,250,249,253,255,255,255,249,249,252,255,230,154,92,84,118,241,251,255;255,255,255,252,249,250,249,249,255,255,255,255,249,249,252,255,255,222,156,144,193,234,248,255;255,255,255,254,249,249,249,251,255,255,255,255,249,249,252,255,252,239,204,197,226,237,250,255;255,255,255,255,250,249,249,253,255,255,255,255,249,249,252,255,255,251,250,249,250,252,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
        AVG.Params.Icons.zoomIcon(:,:,3) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,251,244,251,255;255,255,255,255,255,255,255,255,254,252,252,254,255,255,255,255,255,255,255,255,237,234,252,255;255,238,237,246,237,235,219,160,128,109,121,150,213,254,245,238,234,237,242,237,234,234,234,247;255,234,234,235,238,196,130,125,136,130,120,113,111,162,229,238,238,234,241,238,234,237,241,252;252,234,235,252,220,121,155,237,247,231,218,179,103,71,185,255,242,234,245,245,234,244,255,255;248,234,241,240,143,141,237,255,242,211,204,219,161,77,109,219,237,234,250,241,234,248,255,255;244,234,245,231,114,177,245,247,220,196,194,189,151,105,87,170,234,234,254,237,234,251,255,255;239,234,250,225,94,194,247,228,201,185,187,179,151,113,71,142,234,235,255,235,234,237,251,255;238,234,254,223,81,180,238,217,191,174,186,202,198,114,70,154,237,237,255,246,237,237,254,255;255,255,255,225,104,160,211,205,190,187,213,249,204,113,87,168,255,255,255,255,255,255,255,255;254,248,248,232,129,107,176,202,166,195,230,203,138,85,109,207,255,254,255,255,255,255,255,255;255,249,248,249,208,101,103,137,144,164,164,115,81,80,89,169,252,255,255,255,255,255,255,255;255,252,248,248,252,180,107,90,91,89,90,89,85,81,85,96,184,255,255,255,255,255,255,255;255,254,248,248,252,254,193,125,102,98,107,129,146,99,87,87,96,183,255,255,255,255,255,255;255,255,249,248,250,255,255,244,239,241,245,244,215,141,93,87,87,93,193,249,248,249,252,255;255,255,251,248,249,255,254,244,237,235,235,232,231,216,142,99,87,83,94,191,249,249,253,255;255,255,253,248,248,254,252,248,247,252,251,250,243,245,229,154,98,80,80,85,197,255,255,255;255,255,255,248,248,252,250,248,251,255,255,255,248,248,250,230,146,92,80,80,96,248,255,255;255,255,255,250,248,251,249,248,252,255,255,255,248,248,252,255,230,154,92,83,117,240,250,255;255,255,255,252,248,249,248,248,255,255,255,255,248,248,252,255,255,222,156,144,192,233,246,255;255,255,255,254,248,248,248,250,255,255,255,255,248,248,252,255,252,238,203,196,225,236,250,255;255,255,255,255,249,248,248,252,255,255,255,255,248,248,252,255,255,251,249,248,249,252,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];

        AVG.Params.Icons.syncIcon =...   
           [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
            NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
            NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ;
            NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ;
            NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ;
            NaN,NaN,0  ,0  ,0  ,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ;
            NaN,0  ,0  ,0  ,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN;
            NaN,0  ,0  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN;
            0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,0  ,0  ,0  ,NaN,NaN;
            0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
            0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN;
            0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN;
            NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN;
            NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
        AVG.Params.Icons.syncIcon = repmat(AVG.Params.Icons.syncIcon, [1 1 3]);

        [a,map]=imread('open_files.jpg');
        g=nan(20,20,3);
        g(:)=a(1:5:end,1:5:end,:);
        g(g>128)=NaN;
        g(g<=128)=0;
        AVG.Params.Icons.load=g;
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end %EOF

%% Remarks
%*Filter object does not deal with NaNs. If they appear in the signal the singal becomes all NaN. This causes problems with filtering at the beginning
%   and ending of the session. Zeros can be inserted instead of NaNs in the object that gets the data, but this will result in twisted traces due to the
%   strong transition on onset.