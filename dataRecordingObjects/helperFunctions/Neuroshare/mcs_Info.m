 function [varargout] = mcs_Info(hfile,varargin)
% mcs_Info is a function to show information about the streams that are in a certain file given
% as stream, type, n, TimeSpan and infotext
% 
% Usage:
%               [varargout] = mcs_Info(hFile, varargin)
%
% Description:
%               Generates information about entities. This information is returned in a structured output(table) via CommandWindow 
%               or as a cell matrix of this information.
%
% Parameters:
%               hFile	        Handle/Indentification number(s) to an open file
%               varargin        Variable length input argument list,
%                               "variable argument in"
%
% Return Values:
%               without return value:  information is shown on the display as table
%               with return value:     cell matrix of this information
%
%-------------------------------------------------------------------------------------------------------------------
                        

[nsresult,FileInfo] = ns_GetFileInfo(hfile);                               % Retrieves file information and entity counts

entities = {};                                                             % empty vektor
k=1;

for j = 1:FileInfo.EntityCount                                             % j takes every value between 1 and "Entity Count", saved in "FileInfo"
    
    [nsresult,entity] = ns_GetEntityInfo(hfile,j);                         % retrieves general entity information and type
    stream = entity.EntityLabel(1:8);                                      % definition of "stream" as 1-8 entry of "Entity Label"
                                                                           % The entities are named as follow: 8 letter: stream name with:
                                                                           % 4 letters name
                                                                           % 4 digits stream number
    
    
    h = [stream,int2str(entity.EntityType)];                               % definition of vector h, int2str converts integer to string with integer format 
    entities(k) = {h};                                                     % definition of entities , no longer empty vector
        
    if entity.EntityType == 1                                              % differentation of entity-types by using an if-condition
        infotext(k)= {'Event  Entity'};                                     % if the condition is fulfilled, "infotext" contains "..."
    
    elseif entity.EntityType == 2  
        infotext(k)={'Analog Entity'};                                     % otherwise "infotext" is empty
    
    elseif entity.EntityType == 3  
        infotext(k)={'Analog Entity'};
    end
              
    k = k + 1; 

end

[f,first_oc,dum]= unique(entities);                                        % unique returns the same values as in entities but with no repetitions

for i = 1:length(f)
    n(i) = sum(ismember(entities,f{i}));                                   % ismember returns an array the same size as entities, containing logical 1 (true)
                                                                           % where the elements of entities are in the set f{i} , and logical 0 (false) elsewhere 
end                                                                        % n counts types of different entities

%if nargout==0                                                             % if-condition proofs number of output arguments,nargout returns the number of output
                                                                           % arguments specified in a call to the currently executing function
    
disp('     stream       type      infotext        n       TimeSpan ')      % output "stream","type","n", "TimeSpan" and infotext is given via CommandWindow

%if nargout==0                                                             % if-condition proofs number of output arguments,nargout returns the number of output
                                                                           % arguments specified in a call to the currently executing function

for i = 1:length(f)                          
    

    disp(sprintf('     %s      %d     %s      %d        %f  ', f{i}(1:8) ,   str2num(f{i}(9:9))  ,   infotext{first_oc(i)},  n(i) , FileInfo.TimeSpan )) 
                                                                           % output is given via CommandWindow using sprintf to format data into string
                                                                           % str2num converts string to number 
%end 

%else                                                                      % otherwise
    
for i = 1:length(f)                                                        % output structured as a matrix m
    
    m{i,1} = f{i}(1:8);                                                    % "stream"
    m{i,2} = num2str(str2num(f{i}(9:9)));                                  % "type"
    m{i,3} = infotext{first_oc(i)};                                        % "infotext"
    m{i,4} = n(i);                                                % "n"
    m{i,5} = num2str(FileInfo.TimeSpan);                                   % "TimeSpan"
                                                                           
    
    
end     
    varargout(1)= {m};                                                     %  first entry of varargout is defined as matrix m  
    
    
    
end
end