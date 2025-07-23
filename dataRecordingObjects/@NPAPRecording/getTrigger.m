function [T_ms,chNumber]=getTrigger(obj,startTime_ms,window_ms)
            
    syncdir = string(obj.recordingDir) + filesep + "sync_events"; %Determine dir. 
    
    file = dir (syncdir);
    filenames = {file.name}; %Get filenames in sync folder
    
    filenames = string(filenames(contains(filenames,'out'))); %Get all file names named out. 

    T_ms = cell(1,length(filenames)); %Determine number of cells in T_ms
    
    chNumber = 1:length(T_ms);

    syncdir = string(obj.recordingDir) + filesep + "sync_events"; %Determine dir.

    for n =1:length(filenames)
    
        T_ms{n} = 1000*importdata(string(syncdir) +filesep+ filenames(n))';
    
    end


end
                 
  