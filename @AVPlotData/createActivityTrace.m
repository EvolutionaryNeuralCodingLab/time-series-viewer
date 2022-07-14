function [obj]=createActivityTrace(obj,hControlPanel,hPlotAxis)

obj.refreshPlot=1;

%define default values
obj.plotParams.verticalShift=NaN;
obj.plotParams.plotChannelNumbers=1;
obj.plotParams.plotAnalogData=0;
obj.plotParams.analogScaling=1;

obj.hPlotControls.plotPropGrid=uix.Grid('Parent', obj.hControlPanel, 'Padding', 5, 'Spacing', 10);

obj.hPlotControls.isOverlap=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackNoOverlapPush,'Style','push', 'String','no overlap','string','remove overlap');
obj.hPlotControls.autoScaleYStd=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackAutoScaleStdYPush, 'Style','push', 'String','Auto scale Y (std)');
obj.hPlotControls.verticalShiftEdit=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackVerticalShiftEdit,'Style','edit', 'String','vertical shift');
obj.hPlotControls.plotChannelNumbersCheckbox=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackPlotChannelNumbersCheckbox, 'Style','checkbox','value',obj.plotParams.plotChannelNumbers,'String','plot channel numbers');
obj.hPlotControls.plotChannelNumbersCheckbox=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackPlotAnalogDataCheckbox, 'Style','checkbox','value',obj.plotParams.plotAnalogData,'String','plot analog data');
obj.hPlotControls.analogScalingEdit=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackAnalogScalingEdit,'Style','edit', 'String','analog scaling (usually 0.001)');
set(obj.hPlotControls.plotPropGrid, 'Widths',-2,'Heights', [30 30 30 30 30 30] );

%callback functions for plot controls
    function CallbackPlotAnalogDataCheckbox(hObj,event)
        obj.plotAnalogChannels=obj.hPlotControls.plotChannelNumbersCheckbox.Value;
        %can not replot here since data has to be first extracted and this happens in time series viewer.
    end
    function CallbackAnalogScalingEdit(hObj,event)
        obj.analogScaling=str2double(obj.hPlotControls.analogScalingEdit.String);
        obj.replot;
    end

    function CallbackNoOverlapPush(hObj,event)
        if obj.nCh>1 || obj.nTrials>1
            M=squeeze(obj.M);
        else
            M=squeeze(obj.M)';
        end
        M = double(M);
        minM=nanmin(M,[],2);
        maxM=nanmax(M,[],2);
        obj.plotParams.shifts=flipud(cumsum([0;maxM(end:-1:2)-minM(end-1:-1:1)]));
        set(obj.hPlotControls.verticalShiftEdit,'string','no overlap');
        obj.plotParams.yl=[obj.plotParams.shifts(end)+minM(end) obj.plotParams.shifts(1)+maxM(1)];
        ylim(obj.hPlotAxis,obj.plotParams.yl);
        obj.replot;
    end
    function CallbackVerticalShiftEdit(hObj,event)
        M=squeeze(obj.M);
        verticalShift=max(str2num(get(hObj,'string')),eps);
        obj.plotParams.shifts=((size(M,1)-1):-1:0)'*verticalShift;
        obj.plotParams.yl=[obj.plotParams.shifts(end)-1*verticalShift obj.plotParams.shifts(1)+1*verticalShift];
        ylim(obj.hPlotAxis,obj.plotParams.yl);
        obj.replot;
    end
    function CallbackAutoScaleStdYPush(hObj,event)
        M=squeeze(obj.M);
        verticalShift=nanstd(double(M(:)))*1;
        obj.plotParams.shifts=((size(M,1)-1):-1:0)'*verticalShift;
        set(obj.hPlotControls.verticalShiftEdit,'string',num2str(verticalShift));
        obj.plotParams.yl=[obj.plotParams.shifts(end)-2*verticalShift obj.plotParams.shifts(1)+2*verticalShift];
        ylim(obj.hPlotAxis,obj.plotParams.yl);
        obj.replot;
    end
    function CallbackPlotChannelNumbersCheckbox(hObj,event)
        obj.plotParams.plotChannelNumbers=get(obj.hPlotControls.plotChannelNumbersCheckbox,'value');
        obj.replot;
    end
end %EOF