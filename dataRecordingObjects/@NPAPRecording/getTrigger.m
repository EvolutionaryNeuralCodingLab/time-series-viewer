            
function [T_ms,chNumber]=getTrigger(obj,startTime_ms,window_ms)
            
    syncdir = string(obj.recordingDir) + filesep + "sync_events"; %Determine dir. 
    
    file = dir (syncdir);
    filenames = {file.name}; %Get filenames in sync folder
    
    filenames = string(filenames(contains(filenames,'out'))); %Get all file names named out. 
    
    T_ms = cell(1,round(length(filenames)/2)); %Determine number of cells in T_ms
    
    chNumber = 1:length(T_ms);

    j =1;

    for n = 1:length(filenames)
    
        %CH1
        ch = filenames(1); %Firsdt file name
        dt = 1000*importdata(string(syncdir) +filesep+ string(ch)); %Load file in miliseconds
        filenames = filenames(2:end); %reduce number of filenames by one.
        chi = filenames(contains(filenames,regexp(ch,'\d*','Match'))); %Find the inverse of the channel by selectin an out that matches ch number.
    
        if ~isempty(chi) %If there is a match then combine both the up and down signal of the channel
            dti = 1000*importdata(string(syncdir) +filesep+ string(chi));
            dtc = sort([dt',dti']);

            if ~exist('startTime_ms','var') && ~exist('window_ms','var')

            startTime_ms = 0;
            window_ms = max(dtc)+1;
            end

            dtc = dtc(dtc>= startTime_ms & dtc <= startTime_ms+window_ms);
            T_ms(j) = num2cell(dtc,[1 2]); %asign it to the first position of T_ms cell array
            j = j+1;

        elseif (sum(contains(filenames,regexp(ch,'\d*','Match'))) <1) && ~contains(ch,"inv") %If there is no number match (no inverse signal) and it is not an inverse signal. 
            
            if ~exist('startTime_ms','var') && ~exist('window_ms','var')

            startTime_ms = 0;
            window_ms = max(dt)+1;
            end
            
            dt = dt(dt>= startTime_ms & dt <= startTime_ms+window_ms);
            T_ms(j) = num2cell(dt,[1 2]);
            j = j+1;
        end     

    end
    
end
             
  