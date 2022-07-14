function  [entities] = mcs_GetEntities(hfile,stream_name,varargin)
% mcs_GetEntities is a function to generate a vector of entity numbers of entities with a given characteristic.
% 
% Usage:
%                [entities] = mcs_GetEntities(hFile,stream_name, varargin)
%
% Description:
%                Generates entities out of  data file referenced by 
%                hFile. This information is returned in the structure entities.
%
% Parameters:
%                hFile	        Handle/Indentification number to an open file.
%                stream_name     e.g. 'elec0001'
%                varargin        Variable length input argument list
%                               "variable argument in", (optional) number
%                                of entity type
%
% Return Values:
%                entities       Array of entity numbers
                        

%-------------------------------------------------------------------------------------------------------------------------------------------




[ns_RESULT,FileInfo]=ns_GetFileInfo(hfile);                                % retrieves file information and entity counts, information is returned in the structure "FileInfo"
entities = [];                                                             % entities = empty vector
k=1;


for j= 1:FileInfo.EntityCount                                              % definition j is EntityCount
    
    [ns_RESULT,EntityInfo] = ns_GetEntityInfo(hfile,j);                    % ns_GetEntityInfo retrieves general entity information and type
      
    stream = EntityInfo.EntityLabel(1:length(stream_name));                % definition of "stream"  
     
    if strcmp(stream, stream_name)== 1                                     % strcmp compares two strings for equality. The strings are considered to be equal
                                                                           % if the size and content of each are the same. The function returns a scalar logical
                                                                           % 1 for equality, or scalar logical 0 for inequality.
        
        if length(varargin)==0 || varargin{1}== EntityInfo.EntityType      % if-condition 
            
            entities(k) = j;                                               % definition entities, no longer an empty vector
            
            k = k + 1;
        end
        
    end
    
    
end



end