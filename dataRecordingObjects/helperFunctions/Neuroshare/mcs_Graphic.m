function mcs_Graphic(hfile, EntityID,firstItem,itemCount)
% mcs_Graphic plots analog and triggered entities  
% 
% Usage:
%               mcs_Graphic(hFile, EntityID, firstItem, ItemCount)
%
% Description:
%               Generates a graphical output 
%
% Parameters:
%               hFile	        Handle/Indentification number to an open file.
%               EntityID        Identification number of the entity in the data file
%               firstItem       start
%               itemCount       number of items
%         
% Return Values:
%                /
                   


locX = EntityID;
locY = EntityID;

k=1;
%----------------------------------------------------------------------------------------------------------------------------------------------------------------
for m = EntityID                                                           % identification number of the entity in the data file saved as m
    
    [ns_RESULT,entity] = ns_GetEntityInfo(hfile,m);                        % retrieves general entity information and type by using "hfile" and EntityID = m
                                                                           % the information is passed in the structure entity
 
    
    switch entity.EntityType                                               % switch can differentiate between 
                                                                           % two different kinds of possible input entities 
        
        case 2 % analog
            [ns_RESULT,analog] = ns_GetAnalogInfo(hfile,m);                % "analog" is a structure to receive the Analog Entity
                                                                           % information
     
            locX(k)= analog.LocationX;                                     % locationX is a part of "analog"
            locY(k)= analog.LocationY;                                     % locationY is a part of "analog"
       
        case 3 % segment
            
            [ns_RESULT,segment] = ns_GetSegmentInfo(hfile,m);              % "segment" is a structure that receives segment
                                                                           % information for the requested Segment Entity   
            
            [ns_RESULT,segment_source] = ns_GetSegmentSourceInfo(hfile,m,1);% segment_source is a structure that receives information about the sources that
                                                                           % generated the segment data, implemented by "hfile", EntityID, SourceID
    
            locX(k)= segment_source.LocationX;                             % locationX is a part of "segment_source"
            locY(k)= segment_source.LocationY;                             % locationY is a part of "segment_source"
    end
    k = k + 1;
    
end
%-----------------------------------------------------------------------------------------------------------------------------
min_step_x = 1;                                                            % definition of min_step_x

if (length(unique(locX)) > 1)                                              % if- condition considers length of unique(locX),
                                                                           % unique(locX) contains all entries of locX without any repetitions  
    
    
    min_step_x = min(diff(unique(locX)));                                  % if the condition is fulfilled, min_step_x is changed 
end

tem1=locX-min(locX);                                                       % definition of tem1 
n_step_x =  max(tem1)/min_step_x;                                          % definition of n_step_x

if (length(unique(locX)) > 1)                                              % if- condition considers length of unique(locX),
                                                                           % unique(locX) contains all entries of locX without any repetitions
    
    min_step_x = min_step_x / max(tem1);                                   % if the condition is fulfilled, min_step_x is changed 
    x_Koord = tem1/max(tem1);                                              % definition of x_Koord
    min_step_x= min_step_x * n_step_x/(n_step_x + 1);                      % change of min_step_x
    x_Koord =  x_Koord * n_step_x/(n_step_x + 1);                          % change of x_Koord
else
    min_step_x = 1;                                                        % otherwise min_step_x keeps the value given by first definition
    x_Koord = 0 * tem1;                                                    % and x_Koord changes to 0*tem1
end
%----------------------------------------------------------------------------------------------------------------------------
min_step_y = 1;                                                            % definition of min_step_y
           
if (length(unique(locY)) > 1)                                              % if- condition considers length of unique(locY),
                                                                           % unique(locY) contains all entries of locY without any repetitions 
    
    min_step_y = min(diff(unique(locY)));                                  % if the condition is fulfilled, min_step_y is changed 
end

tem2=locY-min(locY);                                                       % definition of tem2 
n_step_y =  max(tem2)/min_step_y;                                          % definition of n_step_y 

if (length(unique(locY)) > 1)                                              % if- condition considers length of unique(locY),
                                                                           % unique(locY) contains all entries of locX without any repetitions
    
    min_step_y = min_step_y / max(tem2);                                   % if the condition is fulfilled, min_step_y is changed
    y_Koord = tem2/max(tem2);                                              % definition of y_Koord
    min_step_y= min_step_y * n_step_y/(n_step_y + 1);                      % if the condition is fulfilled, min_step_y is changed 
    y_Koord =  y_Koord * n_step_y/(n_step_y + 1);                          % change of y_Koord 
else
    min_step_y = 1;                                                        % otherwise min_step_y keeps the value given by first definition
    y_Koord = 0 * tem2 ;                                                   % and y_Koord changes to 0*tem2
end
%---------------------------------------------------------------------------------------------------------------------------
k=1;                                                 
figure(1)                                                                  % figure creates figure graphics objects                                  

for j=EntityID                                                             % identification number of the entity in the data file saved as j
   
    [ns_RESULT,entity] = ns_GetEntityInfo(hfile,j);                        % retrieves general entity information and type by using "hfile" and EntityID = j
                                                                           % the information is passed in the structure entity
 

    switch entity.EntityType                                               % switch can differentiate between 
                                                                           % two different kinds of possible input entities 
         
        
        case 2
            [ns_RESULT,analog] = ns_GetAnalogInfo(hfile,j);                % "analog" is a structure to receive the Analog Entity
                                                                           % information
                                                          
            x_plus = min_step_x * 0.1;                                     % definition of x_plus, parameter influences the size of axis cross,
                                                                           % choose between 0.1 (large) and 0.4 (small) 
            y_plus = min_step_y * 0.2;                                     % definition of y_plus, parameter influences the size of axis cross,
                                                                           % choose between 0.1 (large) and 0.4 (small)
            
            
            
            
            
            axes('position',[(x_Koord(k)+ x_plus) , (1-y_Koord(k)-min_step_y+y_plus) , (min_step_x-2*x_plus) , (min_step_y-2*y_plus)],'FontSize',8); 
                                                                           % axes creates an axes graphics object in the current figure
                                                                           % using default property values, four element vektor[left bottom width height]
                                                                           % defines position of axes 
      
            fi = firstItem;                                                % definition fi = StartIndex
            ic = itemCount;                                                % definition ic = IndexCount
      
            while (ic ~= 0)                                                % repeatedly execute statements while condition is true
          
                [ns_RESULT,ContCount,data]=ns_GetAnalogData(hfile,j,fi,ic);% retrieves analog data by index using common parameters 
                                                                           % in addition to fi and ic and returns data, number of continuous data points present
                                                                           % in the data as ContCount
      
                max_data = max(data);                                      % definition max_data 
                min_data = min(data);                                      % definition max_data
      
      
                data = data(1:ContCount);                                  % definition data 
                x = fi/analog.SampleRate:1/analog.SampleRate:(fi+ContCount-1)/analog.SampleRate; % definition x
       
                hold on
                plot(x,data)                                               % determines which data is ploted
                hold on
        
                if ic~= ContCount                                          % if- condition considers number of continuous data points
                    plot ([(fi+ContCount)/analog.SampleRate,(fi+ContCount)/analog.SampleRate],[min_data - 0.1*(max_data-min_data), max_data + 0.1*(max_data- min_data)],'-r')
                                                                           % plotts a red mark when a section of continuous data ends
                end
        
                fi = fi + ContCount;                                       % definition fi = StartIndex
                ic = ic - ContCount;                                       % definition ic = IndexCount
        
            end
      
            
           
            xlabel('Time in [s]','fontsize',7)         % the label appears beneath its respective axis in a two-dimensional plot      
            ylabel(['[ '  analog.Units  ' ]'],'fontweight','bold','fontsize',9)
            title(entity.EntityLabel(26:27),'fontsize',7)% the title is located at the top and in the center of the axes
 %-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------       
          case 3
              [ns_RESULT,segment] = ns_GetSegmentInfo(hfile,j);            % "segment" is a structure that receives segment
                                                                           % information for the requested Segment Entity
                                                                 
              
              [ns_RESULT,segment_source] = ns_GetSegmentSourceInfo(hfile,j,1); % Retrieves information about the sources that
                                                                              % generated the segment data
 
                
              x_plus = min_step_x/10;                                      % definition of x_plus
              y_plus = min_step_y/10;                                      % definition of y_plus
              
              axes('position',[(x_Koord(k)+x_plus) , (1-y_Koord(k)-min_step_y+y_plus) , (min_step_x-2*x_plus) , (min_step_y-2*y_plus)],'FontSize',8);
                                                                           % axes creates an axes graphics object in the current figure
                                                                           % using default property values, four element vektor[left bottom width height] 
                                                                           % defines position of axes
                                                                           % for p = firstItem : firstItem + itemCount-1 
              
               for p = fi : fi + ic-1 
                   
                  if p <= entity.ItemCount                                 % if- condition is fulfilled if p is smaller than "ItemCount"
                        
                      [ns_RESULT,timestamp, data, samplecount, unitID] = ns_GetSegmentData(hfile,j,p); % retrieves segment data by index and returns information about 
                                                                                                      % TimeStamp, Data, SampleCount, UnitID	
                      data = data(1:samplecount);                          % workaround for error in function ns_GetSegmentData: data might be longer than samplecount
                                                                           % with invalid data after samplecount

                      x = 0:1/segment.SampleRate:samplecount/segment.SampleRate-1/segment.SampleRate; % definition x
                    
                 
                      hold on 
                      plot(x,data)                                         % plots data versus x 
                
                    end
                end
                
                xlabel('Time in [s]')                                      % the label appears beneath its respective axis in a two-dimensional plot  
                ylabel(['[' segment.Units ']'])
                title(entity.EntityLabel(26:27))                           % the title is located at the top and in the center of the axes
                
      end
      
        k = k + 1;         
         
end

end