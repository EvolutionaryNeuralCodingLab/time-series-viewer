% [hScaleBar]=addScaleBar(hAxis,varargin)
% Function purpose : Adds scale bar to axis
%
% Function recives :    hAxis - axis handle
%                       varagin - 'Property',value
%
% Function give back :  h - axis handle
%
% Last updated : 05/01/10
function [hScaleBar]=addScaleBar(hAxis,varargin)
scaleBarColor='k'; %scale bar color
fontSize=10; %the size of the font on scale bar
scaleBarTransparancy=0.8; %scale bar transparancy
fontColor=[]; %Scale bar font color
XUnitStr='ms'; %x-axis units
YUnitStr='\muV'; %y-axis units
scaleFac=1; %a factor for increasing scale bar size
scaleFacX=[];%a factor for increasing scale bar size on X axis
scaleFacY=[];%a factor for increasing scale bar size on Y axis
xShift=0; %x scale bar shift in units of axis (1 is a shift of a full axis length)
yShift=0; %y scale bar shift in units of axis (1 is a shift of a full axis length)
xOrigin=[]; % The position of scale bar origin on the x axis
yOrigin=[]; % The position of scale bar origin on the y axis
xLim_real=[]; %the limits of the x axis (in XUnitStr units)
yLim_real=[]; %the limits of the y axis (in YUnitStr units)
textOnScaleBar=1; %if==1 text appears on the scale bar, otherwise text appears outside
displayWarnings=0;
scaleBarAxes='xy'; %'x','y','xy'
equalXTWidths=true;
unitsInScaleBar=true;
transparentScale=true;
scaleBarWidth=[]; %in fraction of axes

%print out default arguments and values if no inputs are given
if nargin==0
    defaultArguments=who;
    for i=1:numel(defaultArguments)
        eval(['defaultArgumentValue=' defaultArguments{i} ';']);
        disp([defaultArguments{i} ' = ' num2str(defaultArgumentValue)]);
    end
    return;
end

%Collects all options
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};'])
end

if isempty(scaleFacX)
    scaleFacX=scaleFac;
end
if isempty(scaleFacY)
    scaleFacY=scaleFac;
end

xl=xlim(hAxis);
yl=ylim(hAxis);
if isempty(xLim_real)
    xLim_real=xl; %the limits of the x axis (in XUnitStr units)
end
if isempty(yLim_real)
    yLim_real=yl; %the limits of the y axis (in YUnitStr units)
end

dX=xl(2)-xl(1); %x axis length
dY=yl(2)-yl(1); %y axis length

if isempty(xOrigin)
    xOrigin=0.98*dX+xl(1); %x origin coordinates
end
if isempty(yOrigin)
    yOrigin=0.01*dY+yl(1); %y origin coordinates
end

tmpUnits=get(hAxis,'Units');
set(hAxis,'Units','centimeters');
hPos=get(hAxis,'OuterPosition');
set(hAxis,'Units',tmpUnits);

if isempty(scaleBarWidth)
    if unitsInScaleBar
        scaleBarWidth=0.02;
    else
        scaleBarWidth=0.01;
    end
end

if equalXTWidths
    xWidth=scaleBarWidth*dX*scaleFacX;
    yWidth=scaleBarWidth*dY*(hPos(3)/hPos(4))*scaleFacY;
    %xWidth=scaleBarWidth*dX*scaleFacX;
    %yWidth=scaleBarWidth*dY*scaleFacY*(hPos(3)/hPos(4));
else
    xWidth=scaleBarWidth*dX*scaleFacX;
    yWidth=scaleBarWidth*dY*scaleFacY;
end

dX_real=xLim_real(2)-xLim_real(1);
tmp=floor(log10(dX_real));
xLength_real=10^(tmp-1)*round((dX_real/10^tmp))*scaleFacX;
xLength=xLength_real/(dX_real/dX);

dY_real=yLim_real(2)-yLim_real(1);
tmp=floor(log10(dY_real));
yLength_real=10^(tmp-1)*round((dY_real/10^tmp))*scaleFacY;
yLength=yLength_real/(dY_real/dY);

yTop=yOrigin+yLength;
xLeft=xOrigin-xLength;

if unitsInScaleBar
    xText=[(xOrigin+xLeft)/2+xShift*dX yOrigin+yWidth/2+yShift*dY];
    yText=[xLeft+xLength-xWidth/2+xShift*dX (yOrigin+yTop)/2+yShift*dY];
    alignX='middle';
    alignY='middle';
    if isempty(fontColor)
        fontColor='w';
    end
else
    xText=[(xOrigin+xLeft)/2+xShift*dX yOrigin+yWidth*3/2+yShift*dY];
    yText=[xLeft+xLength+xWidth/2+xShift*dX (yOrigin+yTop)/2+yShift*dY];
    alignX='bottom';
    alignY='top';
    if isempty(fontColor)
        fontColor='k';
    end
end

edgeC=1;
if strcmp(scaleBarAxes,'xy')
    PX=[xOrigin xOrigin xLeft xLeft xOrigin-xWidth xOrigin-xWidth]+xShift*dX;
    PY=[yTop yOrigin yOrigin yOrigin+yWidth yOrigin+yWidth yTop]+yShift*dY;
    %edgeC=[1 2 3 6 5 4];
    xStr=[num2str(xLength_real) XUnitStr];
    yStr=[num2str(yLength_real) YUnitStr];
elseif strcmp(scaleBarAxes,'y')
    PX=[xOrigin xOrigin xOrigin-xWidth xOrigin-xWidth]+xShift*dX;
    PY=[yTop yOrigin yOrigin yTop]+yShift*dY;
    %edgeC=[1 2 3 4];
    xStr='';
    yStr=[num2str(yLength_real) YUnitStr];
elseif strcmp(scaleBarAxes,'x')
    PX=[xOrigin xLeft xLeft xOrigin]+xShift*dX;
    PY=[yOrigin+yWidth yOrigin+yWidth yOrigin yOrigin]+yShift*dY;
    %edgeC=[1 2 3 4];
    xStr=[num2str(xLength_real) XUnitStr];
    yStr='';
    %figure;plot(PX,PY,'.');text(PX,PY,num2str((1:4)'));
end


if any(PX>xl(2) | PX<xl(1) | PY>yl(2) | PY<yl(1))
    if displayWarnings
        disp('Scale bar is outside of axis! changing transparency to make it visible');
    end
    scaleBarOutsideAxis=1;
    scaleBarTransparancy=1;
else
    scaleBarOutsideAxis=0;
end

if transparentScale
    if scaleBarOutsideAxis
        hScaleBar(1,1)=patch(PX,PY,edgeC,'FaceColor',scaleBarColor,'EdgeColor','none','Parent',hAxis,'Clipping','off');
        xlim(xl); %for cases
        ylim(yl);
    else
        hScaleBar(1,1)=patch(PX,PY,edgeC,'FaceColor',scaleBarColor,'EdgeColor','none','FaceAlpha',scaleBarTransparancy,'Parent',hAxis);
    end
else
    lineWidth=4;
    fontColor='k';
    if scaleBarOutsideAxis
        by=line(PX([3 2]),PY([3 2]),'lineWidth',lineWidth,'color',scaleBarColor,'Clipping','off');
        bx=line(PX([3 4]),PY([3 4]),'lineWidth',lineWidth,'color',scaleBarColor,'Clipping','off');
    else
        by=line(PX([3 2]),PY([3 2]),'lineWidth',lineWidth,'color',scaleBarColor);
        bx=line(PX([3 4]),PY([1 2]),'lineWidth',lineWidth,'color',scaleBarColor);
        %error('non transparent scale bar not implemented yet');
    end
    if strcmp(scaleBarAxes,'y')
        delete(bx);
    elseif strcmp(scaleBarAxes,'x')
        delete(by);
    end
end

hScaleBar(2,1)=text(xText(1),xText(2),xStr,...
    'Parent',hAxis,'VerticalAlignment',alignX,'HorizontalAlignment','Center','FontWeight','bold','color',fontColor,'FontSize',fontSize);
hScaleBar(3,1)=text(yText(1),yText(2),yStr,...
    'Parent',hAxis,'Rotation',90,'VerticalAlignment',alignY,'HorizontalAlignment','Center','FontWeight','bold','color',fontColor,'FontSize',fontSize);



