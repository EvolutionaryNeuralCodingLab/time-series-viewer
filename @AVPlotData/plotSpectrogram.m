function [obj]=plotSpectrogram(obj)
%check that input data is valid
if obj.nCh>1 && obj.nTrials>1 
    obj.hPlot=[];hText=[];
    msgbox('spectrogram plot only support either multiple trials or channels (not both)','Attention','error','replace');
    return;
end

Fs=1000/(obj.T(2)-obj.T(1)); %sampling frequency in [Hz]

nRow=ceil(sqrt(obj.nCh*obj.nTrials));
nCol=ceil(obj.nCh*obj.nTrials/nRow);
P=cell(nRow,nCol);
%selection of input data
if obj.nCh==1 && obj.nTrials>1
    for i=1:obj.nCh
        [~,F,T,Ptmp]=spectrogram(squeeze(obj.M(1,i,:)),obj.plotParams.window*Fs/1000,obj.plotParams.overlap*Fs/1000,obj.plotParams.NFFT,Fs);
        P{i}=Ptmp(1:obj.plotParams.maxFreq,:);
    end
elseif obj.nCh>1 && obj.nTrials==1
    for i=1:obj.nCh
        [~,F,T,Ptmp]=spectrogram(squeeze(obj.M(i,1,:)),obj.plotParams.window*Fs/1000,obj.plotParams.overlap*Fs/1000,obj.plotParams.NFFT,Fs);
        P{i}=Ptmp(1:obj.plotParams.maxFreq,:);
    end
elseif obj.nCh==1 && obj.nTrials==1
    i=1;
    [~,F,T,Ptmp]=spectrogram(squeeze(obj.M(1,1,:)),round(obj.plotParams.window*Fs/1000),round(obj.plotParams.overlap*Fs/1000),obj.plotParams.NFFT,Fs);
    P{i}=Ptmp(1:obj.plotParams.maxFreq,:);
end
%initialize combined data matrix
[nFreq nTimes]=size(P{1});
P(i+1:end)={nan([nFreq nTimes])};

M=cell2mat(P);
dT=(T(2)-T(1))*1000;
dF=F(2)-F(1);
%obj.hPlot=surf(repmat(T,[1 nCol]),repmat(F,[nRow 1])',10*log10(abs(M)+eps),'EdgeColor','none');view(0,90);
%obj.hPlot=surf(10*log10(abs(M)+eps),'EdgeColor','none');view(0,90);
obj.hPlot=imagesc(dT/2:dT:(dT*nCol*nTimes),dF/2:dF:(dF*nRow*nFreq),10*log10(M+eps),'Parent',obj.hPlotAxis);

[X,Y]=meshgrid(1:nRow,1:nCol);

if obj.nTrials==1
    hText=text((X(1:obj.nCh)*nTimes-nTimes+1)*dT,(Y(1:obj.nCh)*nFreq-nFreq/8)*dF,num2cell(obj.channelNumbers),...
        'Parent',obj.hPlotAxis,'FontSize',6,'FontWeight','Bold');
elseif obj.nCh==1
    hText=text((X(1:obj.nCh)*nTimes-nTimes+1)*dT,(Y(1:obj.nCh)*nFreq-nFreq/8)*dF,num2cell(1:obj.nTrials),...
        'Parent',obj.hPlotAxis,'FontSize',6,'FontWeight','Bold');
end

hLines=line([([X(1,1:end-1);X(1,1:end-1)])*dT*nTimes [zeros(1,nCol-1);ones(1,nCol-1)*dT*nCol*nTimes]],...
    [[zeros(1,nRow-1);ones(1,nRow-1)*dF*nRow*nFreq] ([Y(1:end-1,1) Y(1:end-1,1)])'*dF*nFreq],...
    'color','k','Parent',obj.hPlotAxis);

xlim(obj.hPlotAxis,[0 (dT*nCol*nTimes)]);
ylim(obj.hPlotAxis,[0 (dF*nRow*nFreq)]);
xlabel(obj.hPlotAxis,'Time [ms]');
ylabel(obj.hPlotAxis,'Frequency [Hz]');
 
set(obj.hPlotControls.spectrogramData,'string',{['T=' num2str(dT,5) ' - ' num2str(T(end)*1000,5)],['F=' num2str(dF,5) ' - ' num2str(F(nFreq),5)]},'FontSize',8);

[hScaleBar]=addScaleBar(obj.hPlotAxis,'YUnitStr','Freq','YUnitStr','Hz');
obj.hPlot=[obj.hPlot;hText;hLines;hScaleBar];