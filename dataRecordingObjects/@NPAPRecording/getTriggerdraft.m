            
function [T_ms,chNumber]=getTrigger(obj,startTime_ms,window_ms)
            
    syncdir = string(NP1.recordingDir) + "\sync_events"; %Determine dir. 
    
    file = dir (syncdir);
    filenames = {file.name}; %Get filenames in sync folder
    
    filenames = string(filenames(contains(filenames,'out'))); %Get all file names named out. 
    
    T_ms = cell(1,round(length(filenames)/2)); %Determine number of cells in T_ms
    
    j =1;


    for n = 1:length(filenames)
    
        %CH1
        ch = filenames(1); %Firsdt file name
        dt = 1000*importdata(string(syncdir) +"\"+ string(ch)); %Load file in miliseconds
        filenames = filenames(2:end); %reduce number of filenames by one.
        chi = filenames(contains(filenames,regexp(ch,'\d*','Match'))); %Find the inverse of the channel by selectin an out that matches ch number.
    
        if ~isempty(chi) %If there is a match then combine both the up and down signal of the channel
            dti = 1000*importdata(string(syncdir) +"\"+ string(chi));
            dtc = sort([dt',dti']);
            T_ms(j) = num2cell(dtc,[1 2]); %asign it to the first position of T_ms cell array
            j = j+1;
        elseif (sum(contains(filenames,regexp(ch,'\d*','Match'))) <1) && ~contains(ch,"inv") %If there is no number match (no inverse signal) and it is not an inverse signal. 
            T_ms(j) = num2cell(dt,[1 2]);
            j = j+1;
        end
    end
    
    if there is a startime and window then: 
    else take whole recording.
    end
    %%Put in cell array
    end
             
                
                
% 
% 
%             if chn == 
%             dti = 
%             ch1i = 'out_1inv';
%             dt1i = 1000*importdata(string(NP.recordingDir) + "\sync_events\" + string(ch1i) + ".txt");
%             
%             if ~exist('startTime_ms','var') && ~exist('window_ms','var')
%                 startTime_ms = min(dt1)-1;
%                 window_ms = max(dt1i)+1;
%             end
%             dt1c = dt1c(dt1c >= startTime_ms & dt1c <= startTime_ms+window_ms);
%             
%             %CH2
%             ch2 = 'out_2';
%             dt2 = 1000*importdata(string(NP.recordingDir) + "\sync_events\" + string(ch2) + ".txt");
%             ch2i = 'out_2inv';
%             dt2i = 1000*importdata(string(NP.recordingDir) + "\sync_events\" + string(ch2i) + ".txt");
%             
%             dt2c = sort([dt2',dt2i']);
%             
%             if ~exist('startTime_ms','var') && ~exist('window_ms','var')
%                 startTime_ms = min(dt2)-1;
%                 window_ms = max(dt2i)+1;
%             end
%             dt2c = dt2c(dt2c >= startTime_ms & dt2c <= startTime_ms+window_ms);
%            
% 
%             ch3 = 'out_3';
%             dt3 = 1000*importdata(string(NP.recordingDir) + "\sync_events\" + string(ch3) + ".txt");
%             if exist('startTime_ms','var') && exist('window_ms','var')
%                 dt3 = dt3(dt3 >= startTime_ms & dt3 <= startTime_ms+window_ms);
%             end
% 
%             T_ms = {dt1c,dt2c,dt3};
% 
%             chNumber = [1,2,3];