function [T_ms,chNumber]=getTrigger(obj,startTime_ms,window_ms)
cutTriggers=false;
if nargin==3
    cutTriggers=true;
elseif nargin==2
    error('Window size is missing from the inputs!')
end

syncdir = string(obj.recordingDir) + filesep + "sync_events"; %Determine dir.

file = dir (syncdir);
filenames = {file.name}; %Get filenames in sync folder

filenames = string(filenames(contains(filenames,'out'))); %Get all file names named out.

T_ms = cell(1,length(filenames)); %Determine number of cells in T_ms

chNumber = 1:length(T_ms);

for n = 1:length(filenames)
    %CH1
    ch = filenames(n); %Firsdt file name
    dt = 1000*importdata(string(syncdir) +filesep+ string(ch)); %Load file in miliseconds

    if ~cutTriggers
        startTime_ms = 0;
        window_ms = max(dt)+1;
    end

    dt = dt(dt>= startTime_ms & dt <= startTime_ms+window_ms);
    T_ms{n} = dt(:)'; %asign it to the first position of T_ms cell array
end
                 
  