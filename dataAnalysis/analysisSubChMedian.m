function analysisSubChMedian(AVPlotDataObj)
AVPlotDataObj.M=bsxfun(@minus,AVPlotDataObj.M,median(AVPlotDataObj.M,1));
