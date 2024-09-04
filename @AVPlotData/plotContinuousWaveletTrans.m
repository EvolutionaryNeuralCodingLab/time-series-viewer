function [obj]=plotContinuousWaveletTrans(obj)
%check that input data is valid
if obj.nCh>1 && obj.nTrials>1 
    obj.hPlot=[];hText=[];
    msgbox('CWT plot only support either multiple trials or channels (not both)','Attention','error','replace');
    return;
end

Fs=1000/(obj.T(2)-obj.T(1)); %sampling frequency in [Hz]


nRow=ceil(sqrt(obj.nCh*obj.nTrials));
nCol=ceil(obj.nCh*obj.nTrials/nRow);
P=cell(nRow,nCol);
%selection of input data
if obj.nCh==1 && obj.nTrials>1
    error('WCT with more than 1 trial was not yet implemented');
    for i=1:obj.nCh
        [cfs,f]=cwt(squeeze(obj.M(1,i,:)),'morse',Fs,'FrequencyLimits',obj.plotParams.freqLim,'VoicesPerOctave',obj.plotParams.voicesPerOctave);
        P{i}=cfs(1:obj.plotParams.maxFreq,:);
    end
elseif obj.nCh>1 && obj.nTrials==1
    error('WCT with more than 1 channel was not yet implemented');
    for i=1:obj.nCh
        [cfs,f]=cwt(squeeze(obj.M(i,1,:)),'morse',Fs,'FrequencyLimits',obj.plotParams.freqLim,'VoicesPerOctave',obj.plotParams.voicesPerOctave);
        P{i}=cfs(1:obj.plotParams.maxFreq,:);
    end
elseif obj.nCh==1 && obj.nTrials==1
    i=1;
    [cfs,f]=cwt(squeeze(obj.M(1,1,:)),'morse',Fs,'FrequencyLimits',obj.plotParams.freqLim,'VoicesPerOctave',obj.plotParams.voicesPerOctave);
    P{i}=cfs;
end

%initialize combined data matrix
[nFreq nTimes]=size(P{1});
P(i+1:end)={nan([nFreq nTimes])};
M=abs(cell2mat(P));

if obj.plotParams.norm2baseline
    pPre=find(obj.T<obj.plotParams.baselineMs,1,'last');
    M=bsxfun(@rdivide,M,mean(M(:,pPre),2));
end

obj.hPlot=imagesc(obj.T,f,M,'Parent',obj.hPlotAxis);
axis(obj.hPlotAxis,'tight'); shading(obj.hPlotAxis,'flat');
xlabel(obj.hPlotAxis,'Time [ms]');
ylabel(obj.hPlotAxis,'Frequency [Hz]');
set(obj.hPlotAxis,'YTick',obj.plotParams.freqLim(1):(obj.plotParams.freqLim(2)-obj.plotParams.freqLim(1))/10:obj.plotParams.freqLim(2),'yscale','log');

%[hScaleBar]=addScaleBar(obj.hPlotAxis,'YUnitStr','Freq','YUnitStr','Hz');
%obj.hPlot=[obj.hPlot;hScaleBar];