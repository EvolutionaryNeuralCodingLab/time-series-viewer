function analysisChannelMedian(AVPlotDataObj)
AVPlotDataObj.M=median(AVPlotDataObj.M,1);
AVPlotDataObj.channelNumbers=0;
AVPlotDataObj.channelNames={'Ch. Med.'};