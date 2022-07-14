function analysisHilberEnvelop(AVPlotDataObj)
if AVPlotDataObj.nTrials==1
    AVPlotDataObj.M(:,1,:)=abs(hilbert(squeeze(AVPlotDataObj.M)'))';
elseif AVPlotDataObj.nCh==1
    AVPlotDataObj.M(1,:,:)=abs(hilbert(squeeze(AVPlotDataObj.M)'))';
else
    for i=1:AVPlotDataObj.nTrials
        AVPlotDataObj.M(:,i,:)=abs(hilbert(squeeze(AVPlotDataObj.M(:,i,:))'))';
    end
end
mean(AVPlotDataObj.M,1);
