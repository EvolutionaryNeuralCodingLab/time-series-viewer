function analysisRemoveChOffset(AVPlotDataObj)
AVPlotDataObj.M=bsxfun(@minus,AVPlotDataObj.M,median(AVPlotDataObj.M(:,:,1:100),3));