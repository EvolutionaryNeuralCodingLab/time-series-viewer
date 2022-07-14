function obj=createSuperimposedActivityTrace(obj,hControlPanel,hPlotAxis)

obj.refreshPlot=1;

%create the GUI plot controls
obj.hPlotControls.plotPropGrid=uix.Grid('Parent', obj.hControlPanel, 'Padding', 5, 'Spacing', 10);
obj.hPlotControls.autoScaleY=uicontrol('Parent', obj.hPlotControls.plotPropGrid, 'Callback',{@CallbackAutoScaleYPush}, 'Style','push', 'String','Auto scale Y');
obj.hPlotControls.autoScaleYStd=uicontrol('Parent', obj.hPlotControls.plotPropGrid, 'Callback',{@CallbackAutoScaleStdYPush}, 'Style','push', 'String','Auto scale Y (std)');

set(obj.hPlotControls.plotPropGrid, 'Widths',-1,'Heights', [30 30] );

%callback functions for plot controls
    function CallbackAutoScaleYPush(hObj,event)
        minVal=min(obj.M(:));
        maxVal=max(obj.M(:));
        range=maxVal-minVal;
        ylim(obj.hPlotAxis,[minVal-0.01*range maxVal+0.01*range]);
    end
    function CallbackAutoScaleStdYPush(hObj,event)
        meanM=mean(obj.M(:));
        stdM=std(obj.M(:));
        ylim(obj.hPlotAxis,[meanM-5*stdM meanM+5*stdM]);
    end

end %EOF