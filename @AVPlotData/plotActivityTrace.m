function [obj]=plotActivityTrace(obj)
%selection of input data
hText=[];hTrigs=[];
if obj.nCh==1 && obj.nTrials>1
    M=squeeze(obj.M);
    %this option fails only when the two dimensions of M are equal - check
elseif obj.nTrials==1 && obj.nCh>1
    M=squeeze(obj.M);
elseif obj.nCh==1 && obj.nTrials==1
    M=squeeze(obj.M)';
elseif obj.nCh>1 && obj.nTrials>1
    obj.hPlot=[];
    msgbox('activity trace plot only support either multiple trials or channels (not both)','Attention','error','replace');
    return;
end

if obj.refreshPlot %for first appearance when shifts were not chosen
    verticalShift=max(nanstd(double(obj.M(:)))*5,eps);
    obj.plotParams.shifts=((size(M,1)-1):-1:0)'*verticalShift;
    set(obj.hPlotControls.verticalShiftEdit,'string',num2str(verticalShift));
    obj.plotParams.yl=[obj.plotParams.shifts(end)-1*verticalShift obj.plotParams.shifts(1)+1*verticalShift+eps];
    ylim(obj.hPlotAxis,obj.plotParams.yl);
    obj.refreshPlot=0;
end

M=bsxfun(@plus,double(M),obj.plotParams.shifts);
obj.hPlot=plot(obj.hPlotAxis,obj.T,M);obj.hPlotAxis.ColorOrderIndex=1;

if obj.nTrials==1 & obj.plotParams.plotChannelNumbers
    hText=text(obj.T(end)*0.99*ones(1,obj.nCh),M(:,end),num2cell(obj.channelNumbers),...
        'Parent',obj.hPlotAxis,'FontSize',6,'FontWeight','Bold','BackgroundColor','w');
elseif obj.nCh==1 & obj.plotParams.plotChannelNumbers
    hText=text(obj.T(end)*0.99*ones(1,obj.nTrials),M(:,end),num2cell(1:obj.nTrials),...
        'Parent',obj.hPlotAxis,'FontSize',6,'FontWeight','Bold','BackgroundColor','w');
end

%analog channels are added as regular electrode channels (maybe worth changing in the future).
nA=size(obj.A,1);
if obj.plotAnalogChannels && nA>0
    A=squeeze(obj.A);
    if nA==1
        A=A';
    end
    yl=obj.plotParams.yl;
    dV=(yl(2)-yl(1))/obj.nCh;
    A=bsxfun(@plus,bsxfun(@minus,double(A),mean(A,2))*obj.analogScaling,(yl(1)-dV*(1:nA))');
    
    obj.hPlot=[obj.hPlot;plot(obj.hPlotAxis,obj.TA,A,'k')];
    hText=[hText;text(obj.T(end)*0.98*ones(1,nA),A(:,end),obj.analogChannelNames,...
        'Parent',obj.hPlotAxis,'FontSize',6,'FontWeight','Bold','BackgroundColor','w')];
    
    ylim(obj.hPlotAxis,[yl(1)-dV*(nA+1) yl(2)]);
else
    ylim(obj.hPlotAxis,obj.plotParams.yl);
end

if ~isempty(obj.trigMarks)
    hTrigs=line(obj.hPlotAxis,[obj.trigMarks;obj.trigMarks]',obj.hPlotAxis.YLim','Color',[0.8 0.8 0.8]);
end

xlim(obj.hPlotAxis,[obj.T(1) obj.T(end)]);
xlabel(obj.hPlotAxis,'Time [ms]','FontSize',14);
ylabel(obj.hPlotAxis,'Voltage [\muV]','FontSize',14);
[hScaleBar]=addScaleBar(obj.hPlotAxis);
obj.hPlot=[obj.hPlot;hText;hScaleBar;hTrigs];