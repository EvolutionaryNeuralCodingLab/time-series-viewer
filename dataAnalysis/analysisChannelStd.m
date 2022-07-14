function analysisChannelStd(AVPlotDataObj)
AVPlotDataObj.M=std(AVPlotDataObj.M,1);
AVPlotDataObj.channelNumbers=0;
AVPlotDataObj.channelNames={'Ch. Std.'};