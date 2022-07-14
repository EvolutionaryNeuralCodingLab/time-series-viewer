function [obj]=plotRasterPlot(obj)
hText=[];
if obj.refreshPlot
    obj.refreshPlot=0;
end
Ydev=obj.nCh;
ylab='Channel';
if obj.nCh==1 && obj.nTrials>1
    M=squeeze(obj.M);
    Ydev=obj.nTrials;
    ylab='Trial';
elseif obj.nTrials==1 && obj.nCh>1
    M=squeeze(obj.M);
elseif obj.nCh==1 && obj.nTrials==1
    M=squeeze(obj.M)';
elseif obj.nCh>1 && obj.nTrials>1
    obj.hPlot=[];hText=[];
    msgbox('activity trace plot only support either multiple trials or channels (not both)','Attention','error','replace');
    return;
end

if obj.plotParams.constantBinNumber %smear data on more bins to make spikes wider (so that they are always seen)
    convRatio=ceil(obj.nSamples/obj.plotParams.binNumber);
    M=convn(squeeze(M),ones(1,convRatio),'same');
    M=M(:,1:convRatio:end);
end
imageMax=max(M,[],2);
if obj.plotParams.normManual
    M=M/obj.plotParams.maxValue;
else
    if obj.plotParams.normChannel %normalize each channel to its local max
        for i=1:size(M,1)
            M(i,:)=M(i,:)./imageMax(i);
        end
    else %normalize to one global maximum
        M=M/max(imageMax);
    end
end

%plot electrode raster
obj.hPlot=image(obj.T,1:Ydev,M,'CDataMapping','scaled','Parent',obj.hPlotAxis);
set(obj.hPlotAxis,'CLim',[0 1],'YDir','reverse');
  
xlim(obj.hPlotAxis,[0.5 obj.T(end)+0.5]);
ylim(obj.hPlotAxis,[0.5 Ydev+0.5]);
set(obj.hPlotAxis,'Box','on');

xlabel(obj.hPlotAxis,'Time [ms]');
ylabel(obj.hPlotAxis,ylab);

if obj.nTrials==1 & obj.plotParams.plotChannelNumbers
    hText=text(obj.T(end)*0.99*ones(1,obj.nCh),1:obj.nCh,obj.channelNames,...
    'Parent',obj.hPlotAxis,'FontSize',6,'FontWeight','Bold','BackgroundColor','w');
elseif obj.nCh==1 & obj.plotParams.plotChannelNumbers
    hText=text(obj.T(end)*0.99,round(obj.nTrials/2),obj.channelNames,...
    'Parent',obj.hPlotAxis,'FontSize',6,'FontWeight','Bold','BackgroundColor','w');
end

obj.hPlot=[obj.hPlot;hText];
end
