function [obj]=plotActivityTracePhysicalSpace(obj)
if obj.refreshPlot
    obj.plotParams.minV=min(min(obj.M(:)),-eps);
    obj.plotParams.maxV=max(max(obj.M(:)),eps);
    set(obj.hPlotControls.voltageScaleEdit,'string',num2str(obj.plotParams.maxV-obj.plotParams.minV));
    obj.refreshPlot=0;
end
%Mtmp=(squeeze(obj.M)'-obj.plotParams.minV)/(obj.plotParams.maxV-obj.plotParams.minV)-0.5; %scaling
Mtmp=bsxfun(@rdivide,  bsxfun(@minus,squeeze(obj.M),obj.plotParams.minV)  ,  (obj.plotParams.maxV-obj.plotParams.minV))' - 0.5; %scaling
if size(Mtmp,1)==1
    Mtmp=Mtmp'; %fix this to be independent of size
end

WaveformShiftX=((1/obj.nSamples):(1/obj.nSamples):1)';
obj.hPlot=plot(ones(obj.nSamples,1)*obj.plotParams.x(obj.channelNumbers)'-0.5+WaveformShiftX*ones(1,obj.nCh),...
    ones(obj.nSamples,1)*obj.plotParams.y(obj.channelNumbers)'+Mtmp,'Parent',obj.hPlotAxis,'color',[0.1 0.1 0.9]);

[yE,xE]=size(obj.plotParams.En);
[hScaleBar]=addScaleBar(obj.hPlotAxis,...
    'xLim_real',xE*obj.T([1 end]),'yLim_real',[0 yE*(obj.plotParams.maxV(1)-obj.plotParams.minV(1))]);
obj.hPlot=[obj.hPlot;hScaleBar];
end
