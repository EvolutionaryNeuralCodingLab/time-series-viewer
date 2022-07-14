function [obj]=createActivityTracePhysicalSpace(obj)
%electrode layouts, notice that the first layout of a given number of channels is the default
%load all layouts
electrodeFolder=dir([fileparts(which('layout_200_12x12.mat')) filesep 'layout_*.mat']);
obj.plotParams.allElecLayoutNames=cellfun(@(x) x(1:end-4),{electrodeFolder.name},'UniformOutput',0);
obj.plotParams.elecLayoutNames4Popup=cellfun(@(x) x(8:end),obj.plotParams.allElecLayoutNames,'UniformOutput',0);

if ~isempty(obj.recordingObjReference)
    pLayout=find(strcmp(['layout_',obj.recordingObjReference.layoutName],obj.plotParams.allElecLayoutNames));
else
    pLayout=[];
end
if ~isempty(pLayout)
    obj.plotParams.electrodeLayoutNumber=pLayout;
else
    if ~isempty(obj.nCh)
        elecNums=cellfun(@(x) x(2),regexp(obj.plotParams.allElecLayoutNames,'_','split'));
        elecNums=str2double(elecNums);
        obj.plotParams.electrodeLayoutNumber=find(elecNums>=obj.nCh,1,'first');
    else
        obj.plotParams.electrodeLayoutNumber=1;
    end
end

obj.plotParams.channelDisplayOptions={'shift 2 ch names','shift 2 ch numbers'};
obj.plotParams.channelDisplayCurrentOption=1;
%define default values
obj.refreshPlot=1;

%create the GUI plot controls
obj.hPlotControls.plotPropGrid=uix.Grid('Parent', obj.hControlPanel, 'Padding', 5, 'Spacing', 10);

obj.hPlotControls.noOverlap=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackNoOverlapPush,'Style','push', 'String','no overlap','string','remove overlap');
obj.hPlotControls.autoScaleYStd=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackAutoScaleStdYPush, 'Style','push', 'String','Auto scale Y (std)');
obj.hPlotControls.normalizedScale=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackNormalizedScalePush, 'Style','push', 'String','Normalize scale');
obj.hPlotControls.voltageScaleEdit=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackVerticalShiftEdit,'Style','edit', 'String','vertical shift');
obj.hPlotControls.channelNumberNameSelectionPush=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackChannelNumberNameSelectionPush, 'Style','push', 'String',obj.plotParams.channelDisplayOptions{obj.plotParams.channelDisplayCurrentOption});
obj.hPlotControls.drawGridPush=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackDrawGridCheck,'Style','checkbox', 'String','Draw grid','value',1);
obj.hPlotControls.drawNumbersPush=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackDrawNumbersCheck,'Style','checkbox', 'String','Draw numbers','value',1);
obj.hPlotControls.changeElectrodeLayout=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackChangeElectrodeLayoutPop,'Style','popup', 'String',obj.plotParams.elecLayoutNames4Popup,...
    'value',obj.plotParams.electrodeLayoutNumber);

hold(obj.hPlotAxis,'on');
set(obj.hPlotControls.plotPropGrid, 'Widths',-1,'Heights', [25 25 25 25 25 25] );

obj=preparePhysicalElectrodeSpace(obj);

%nested functions
    function obj=preparePhysicalElectrodeSpace(obj)
        %Scaling parameter in case the data should be super imposed on a figure. see ConnectivityPhysicalSpacePlot
        obj.plotParams.aX=1;obj.plotParams.aY=1;obj.plotParams.bX=0;obj.plotParams.bY=0;
        
        %creating translation between electrodes and locations
        load(obj.plotParams.allElecLayoutNames{obj.plotParams.electrodeLayoutNumber},'En','Ena');
        %for nicer ploting remove rows or columns with all nans,
        p1=all(isnan(En));
        p2=all(isnan(En'));
        En(:,p1)=[];En(p2,:)=[];Ena(:,p1)=[];Ena(p2,:)=[];
        pNotNans=find(~isnan(En));
        [obj.plotParams.allChannels p]=sort(En(pNotNans));
        obj.plotParams.allChannelNames=Ena(pNotNans(p));
        obj.plotParams.nAllChannels=numel(obj.plotParams.allChannels); %changed 2.8.16 - Mark to enable gaps in channel numbers in layout
        %obj.plotParams.nAllChannels=max(obj.plotParams.allChannels);

        if numel(obj.channelNumbers)~=obj.plotParams.nAllChannels
            %msgbox('The number of data channels does not fit layout channels!');
            answer = questdlg('Should I try to ignore channels from the recording?','Electrode missmatch issue','yes','no','no');
            switch answer
                case 'yes'
                case 'no'
                    return;
            end
        end
        
        obj.plotParams.translation=NaN(obj.plotParams.nAllChannels,3);
        for i=1:obj.plotParams.nAllChannels
            [n,m]=find(En==obj.plotParams.allChannels(i));
            obj.plotParams.translation(i,:)=[m n i];
        end
        [max_grid_y,max_grid_x]=size(En);
        
        XL=[0.5 max_grid_x+0.5];
        YL=[0.5 max_grid_y+0.5];
        
        xlim(obj.hPlotAxis,XL);
        ylim(obj.hPlotAxis,YL);
                
        xM=(1.5:(max_grid_y-0.5))';
        yM=(1.5:(max_grid_x-0.5))';
        
        if ~isempty(xM)
            hGrid1=line([XL(1) XL(2)],[xM*obj.plotParams.aY+obj.plotParams.bY xM*obj.plotParams.aY+obj.plotParams.bY],...
                'LineWidth',2,'Color',[0.9 0.9 0.9],'Parent',obj.hPlotAxis);
        else
            hGrid1=[];
        end
        if ~isempty(yM)
            hGrid2=line([yM *obj.plotParams.aX+obj.plotParams.bX yM*obj.plotParams.aX+obj.plotParams.bX],[YL(1) YL(2)],...
                'LineWidth',2,'Color',[0.9 0.9 0.9],'Parent',obj.hPlotAxis);
        else
            hGrid2=[];
        end
        obj.hPlotBackground.hGrid=[hGrid1;hGrid2];
        
        obj.plotParams.x=obj.plotParams.aX*obj.plotParams.translation(:,1)+obj.plotParams.bX;
        obj.plotParams.y=obj.plotParams.aY*obj.plotParams.translation(:,2)+obj.plotParams.bY;
        
        set(obj.hPlotAxis,'xtick',[]);
        set(obj.hPlotAxis,'ytick',[]);
        set(obj.hPlotAxis,'XColor',[.8 .8 .8]); %To avoid axis being seen when printing
        set(obj.hPlotAxis,'YColor',[.8 .8 .8]); %To avoid axis being seen when printing
        
        %Plotting numbers on top of the propagation points
        obj.hPlotBackground.hElectrodes=text(obj.plotParams.aX*(obj.plotParams.translation(:,1)-0.44)+obj.plotParams.bX,0.35...
            +obj.plotParams.aY*(obj.plotParams.translation(:,2))+obj.plotParams.bY,num2cell(obj.plotParams.allChannels),...
            'fontsize',8,'FontWeight','Bold','Parent',obj.hPlotAxis);
        
        obj.plotParams.En=En;
    end

%callback functions for plot controls
    function CallbackChannelNumberNameSelectionPush(hObj,event)
        delete(obj.hPlotBackground.hElectrodes);
        if obj.plotParams.channelDisplayCurrentOption==1
            obj.hPlotBackground.hElectrodes=text(obj.plotParams.aX*(obj.plotParams.translation(:,1)-0.44)+obj.plotParams.bX,0.35...
                +obj.plotParams.aY*(obj.plotParams.translation(:,2))+obj.plotParams.bY,obj.plotParams.allChannelNames,...
                'fontsize',8,'FontWeight','Bold','Parent',obj.hPlotAxis);
            obj.plotParams.channelDisplayCurrentOption=2;
        else
            obj.hPlotBackground.hElectrodes=text(obj.plotParams.aX*(obj.plotParams.translation(:,1)-0.44)+obj.plotParams.bX,0.35...
                +obj.plotParams.aY*(obj.plotParams.translation(:,2))+obj.plotParams.bY,num2cell(obj.plotParams.allChannels),...
                'fontsize',8,'FontWeight','Bold','Parent',obj.hPlotAxis);
            obj.plotParams.channelDisplayCurrentOption=1;
        end
        set(obj.hPlotControls.channelNumberNameSelectionPush,'string',obj.plotParams.channelDisplayOptions{obj.plotParams.channelDisplayCurrentOption});
    end
    function CallbackChangeElectrodeLayoutPop(hObj,event)
        obj.plotParams.electrodeLayoutNumber=get(hObj,'Value');
        delete(obj.hPlotBackground.hElectrodes);
        obj=preparePhysicalElectrodeSpace(obj);
        obj.replot;
    end
    function CallbackDrawGridCheck(hObj,event)
        if get(hObj,'value')
            set(obj.hPlotBackground.hGrid,'visible','on');
        else
            set(obj.hPlotBackground.hGrid,'visible','off');
        end
    end
    function CallbackDrawNumbersCheck(hObj,event)
        if get(hObj,'value')
            set(obj.hPlotBackground.hElectrodes,'visible','on');
        else
            set(obj.hPlotBackground.hElectrodes,'visible','off');
        end
    end
    function CallbackNoOverlapPush(hObj,event)
        obj.plotParams.minV=min(obj.M(:));
        obj.plotParams.maxV=max(obj.M(:));
        obj.replot;
        set(obj.hPlotControls.voltageScaleEdit,'string',num2str(obj.plotParams.maxV-obj.plotParams.minV));
    end
    function CallbackVerticalShiftEdit(hObj,event)
        voltageScale=str2num(get(obj.hPlotControls.voltageScaleEdit,'string'));
        obj.plotParams.minV=-voltageScale/2;
        obj.plotParams.maxV=voltageScale/2;
        obj.replot;
    end
    function CallbackAutoScaleStdYPush(hObj,event)
        stdV=std(obj.M(:));
        meanV=mean(obj.M(:));
        obj.plotParams.minV=meanV-10*stdV;
        obj.plotParams.maxV=meanV+10*stdV;
        set(obj.hPlotControls.voltageScaleEdit,'string',num2str(obj.plotParams.maxV-obj.plotParams.minV));
        obj.replot;
    end
    function CallbackNormalizedScalePush(hObj,event)
        obj.plotParams.maxV=max(obj.M,[],3);
        obj.plotParams.minV=min(obj.M,[],3);
        set(obj.hPlotControls.voltageScaleEdit,'string','normalized scale');
        obj.replot;
    end

end %EOF