function [obj]=createSpectrogram(obj,hControlPanel,hPlotAxis)

%define default values

%create the GUI plot controls
obj.hPlotControls.plotPropGrid=uix.Grid('Parent', obj.hControlPanel, 'Padding', 10, 'Spacing', 10);

obj.hPlotControls.voicesPerOctaveTxt=uicontrol('Parent', obj.hPlotControls.plotPropGrid,'Style','text', 'String','Voices/Octave','HorizontalAlignment','left');
obj.hPlotControls.baselineMsTxt=uicontrol('Parent', obj.hPlotControls.plotPropGrid,'Style','text', 'String','baseline [ms]','HorizontalAlignment','left');
obj.hPlotControls.norm2baselineTxt=uicontrol('Parent', obj.hPlotControls.plotPropGrid,'Style','text', 'String','Normalize','HorizontalAlignment','left');
obj.hPlotControls.freqLimitsTxt=uicontrol('Parent', obj.hPlotControls.plotPropGrid,'Style','text', 'String','Freq. limits [x,y]','HorizontalAlignment','left');
obj.hPlotControls.tmpTxt=uicontrol('Parent', obj.hPlotControls.plotPropGrid,'Style','text', 'String','','HorizontalAlignment','left');

obj.hPlotControls.voicesPerOctaveEdit=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackVoicesPerOctaveEdit,'Style','edit', 'String','10');
obj.hPlotControls.baselineMsEdit=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackBaselineMsEdit,'Style','edit', 'String','500');
obj.hPlotControls.norm2baselineCheck=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackNorm2baselineCheck,'Style','checkbox', 'Value',true);
obj.hPlotControls.freqLimEdit=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackFreqLimEdit,'Style','edit', 'String','0 200');
obj.hPlotControls.replot=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackReplotPush,'Style','push', 'String','Replot');

set(obj.hPlotControls.plotPropGrid, 'Widths',[-1 -1],'Heights', [30 30 30 30 30] );

%get initial parameters
CallbackVoicesPerOctaveEdit;
CallbackBaselineMsEdit;
CallbackNorm2baselineCheck;
CallbackFreqLimEdit;

%callback functions for plot controls
    function CallbackVoicesPerOctaveEdit(hObj,event)
        obj.plotParams.voicesPerOctave=str2double(get(obj.hPlotControls.voicesPerOctaveEdit,'string'));
    end
    function CallbackBaselineMsEdit(hObj,event)
        obj.plotParams.baselineMs=str2double(get(obj.hPlotControls.baselineMsEdit,'string'));
    end
    function CallbackNorm2baselineCheck(hObj,event)
        obj.plotParams.norm2baseline=obj.hPlotControls.norm2baselineCheck.Value;
    end
    function CallbackFreqLimEdit(hObj,event)
        obj.plotParams.freqLim=round(str2num(obj.hPlotControls.freqLimEdit.String));
        maxAllowedFreq=10000;
        if obj.plotParams.freqLim(2)>maxAllowedFreq
            obj.plotParams.freqLim(2)=maxAllowedFreq;
            msgbox('Max freq too high','Attention','error','replace');
        end
        if obj.plotParams.freqLim(1)<0
            obj.plotParams.freqLim(1)=0;
            msgbox('Min freq can not be smaller than 0','Attention','error','replace');
        end
        obj.hPlotControls.freqLimEdit.String=num2str(obj.plotParams.freqLim);
    end
   
    function CallbackReplotPush(hObj,event)
        obj.replot;
    end
end %EOF