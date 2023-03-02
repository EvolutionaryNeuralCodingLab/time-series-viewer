function analysisStdThreshold(AVPlotDataObj)
AVPlotDataObj.M=std(AVPlotDataObj.M,2);
for i=1:size(AVPlotDataObj.M,2)
    AVPlotDataObj.M=bsxfun(@minus,AVPlotDataObj.M,std(AVPlotDataObj.M,2))
end