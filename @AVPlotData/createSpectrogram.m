function [obj]=createSpectrogram(obj,hControlPanel,hPlotAxis)

%define default values

%create the GUI plot controls
obj.hPlotControls.plotPropGrid=uix.Grid('Parent', obj.hControlPanel, 'Padding', 10, 'Spacing', 10);

obj.hPlotControls.windowTxt=uicontrol('Parent', obj.hPlotControls.plotPropGrid,'Style','text', 'String','Window [ms]','HorizontalAlignment','left');
obj.hPlotControls.overlapTxt=uicontrol('Parent', obj.hPlotControls.plotPropGrid,'Style','text', 'String','Overlap [ms]','HorizontalAlignment','left');
obj.hPlotControls.NFFTTxt=uicontrol('Parent', obj.hPlotControls.plotPropGrid,'Style','text', 'String','NFFT','HorizontalAlignment','left');
obj.hPlotControls.maxFreqTxt=uicontrol('Parent', obj.hPlotControls.plotPropGrid,'Style','text', 'String','max Freq.','HorizontalAlignment','left');
obj.hPlotControls.spectrogramData=uicontrol('Parent', obj.hPlotControls.plotPropGrid,'Style','text', 'String','','HorizontalAlignment','left');

obj.hPlotControls.windowEdit=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackWindowEdit,'Style','edit', 'String','100');
obj.hPlotControls.OverlapEdit=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackOverlapEdit,'Style','edit', 'String','50');
obj.hPlotControls.NFFTEdit=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackNFFTEdit,'Style','edit', 'String','128');
obj.hPlotControls.maxFreqEdit=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackMaxFreqEdit,'Style','edit', 'String','65');
obj.hPlotControls.replot=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackReplotPush,'Style','push', 'String','Replot');


set(obj.hPlotControls.plotPropGrid, 'Widths',[-1 -1],'Heights', [30 30 30 30 30] );

%get initial parameters
CallbackWindowEdit;
CallbackOverlapEdit;
CallbackNFFTEdit;
CallbackMaxFreqEdit;

%callback functions for plot controls
    function CallbackWindowEdit(hObj,event)
        obj.plotParams.window=str2double(get(obj.hPlotControls.windowEdit,'string'));
    end
    function CallbackOverlapEdit(hObj,event)
        obj.plotParams.overlap=str2double(get(obj.hPlotControls.OverlapEdit,'string'));
    end
    function CallbackNFFTEdit(hObj,event)
        obj.plotParams.NFFT=round(str2double(get(obj.hPlotControls.NFFTEdit,'string')));
    end
    function CallbackMaxFreqEdit(hObj,event)
        tmpFreq=round(str2double(get(obj.hPlotControls.maxFreqEdit,'string')));
        maxAllowedFreq=obj.plotParams.NFFT/2+1;
        if tmpFreq<=maxAllowedFreq
            obj.plotParams.maxFreq=tmpFreq;
        else
            obj.plotParams.maxFreq=maxAllowedFreq;
            set(obj.hPlotControls.maxFreqEdit,'string',num2str(maxAllowedFreq));
            msgbox('The maximal frequency number is NFFT/2+1','Attention','error','replace');
        end
    end
    function CallbackReplotPush(hObj,event)
        obj.replot;
    end
end %EOF