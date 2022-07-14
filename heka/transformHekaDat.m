function varargout = transformHekaDat(varargin)

% The original .dat-files are first read using the ImportHEKA function of
% the sigTOOL toolbox (http://sourceforge.net/projects/sigtool/). This functions of
% the toolbox have to be on the Matlab path. Thanks to Malcolm
% Lidierth for help with the sigTOOL functions!
% These data and some of the metadata are then
% transformed into a tree structure and saved as a mat-file for later use.
% ATTENTION: Two errors are in the SigTool functions:
% 1. Function ImportHEKA: Change line 646 to:
%    % Fix added 28.11.2011 - may be 1 char only
%    if numel(tr_s(1).TrXUnit)<2
%        tr_s(1).TrXUnit(2)=' ';
%    end
% 2. Function scCreateKCLFile.m: Change line 43 to: 
%    matfilename=[targetpath filesep name '.kcl'];
%
% Structures created by this function can be visualized using the traceVis
% GUI.

if nargin == 0
    fName = uigetfile('*.dat','Choose HEKA data file');
    if fName == 0
        return
    end
else
    fName = varargin{1};
end

[pn fn ext] = fileparts(fName);
% Load data file
fn = ImportHEKA(fName);
[pn fn ext] = fileparts(fn);
movefile([fn ext], [fn '.mat']);
in = load([fn '.mat']);
% Rename file back to kcl file
movefile([fn '.mat'], [fn '.kcl']);
fName = which(fName);


% Sort data

for i = 1:size(in.FileSource.header,2)
    ind{i} = [];
    for j = 1:size(in.FileSource.header,1)
        a = in.FileSource.header;
        a = a{j,i};
        if ~isempty(a)
            ind{i}(end+1) = j;
        end
    end
end

% Find group indices (e.g. which series belongs to which group; a new group
% can be created by the experimenter by clicking "New group" or "New 
% experiment" in the file menu).
% All series are numbered consecutively starting from 1 and "grouped"
% in a cell array according to the group they belong to
if numel(ind{2}) == 1
    groups{1} = 1:numel(ind{3});
else
    for g = 1:numel(ind{2})
        if g == 1
            groups{g} = 1:sum(ind{3} > ind{2}(g) & ind{3} < ind{2}(g+1));
        elseif g < numel(ind{2})
            groups{g} = groups{g-1}(end) + (1:sum(ind{3} > ind{2}(g) & ind{3} < ind{2}(g+1)));
        else
            groups{g} = groups{g-1}(end) + (1:sum(ind{3} > ind{2}(g)));
        end
    end
end

% Find series indices (e.g. which sweeps belongs to which series)
% All sweeps are numbered consecutively starting from 1 and "grouped"
% in a cell array according to the series they belong to
if numel(ind{3}) == 1
    series{1} = 1:numel(ind{4});
else
    for g = 1:numel(ind{3})
        if g == 1
            series{g} = 1:sum(ind{4} > ind{3}(g) & ind{4} < ind{3}(g+1));
        elseif g < numel(ind{3})
            series{g} = series{g-1}(end) + (1:sum(ind{4} > ind{3}(g) & ind{4} < ind{3}(g+1)));
        else
            series{g} = series{g-1}(end) + (1:sum(ind{4} > ind{3}(g)));
        end
    end
end


% Find sweep indices (e.g. which traces belongs to which sweeps)
% All traces are numbered consecutively starting from 1 and "grouped"
% in a cell array according to the sweep they belong to
if numel(ind{4}) == 1
    sweeps{1} = 1:numel(ind{5});
else
    for g = 1:numel(ind{4})
        if g == 1
            sweeps{g} = 1:sum(ind{5} > ind{4}(g) & ind{5} < ind{4}(g+1));
        elseif g < numel(ind{4})
            sweeps{g} = sweeps{g-1}(end) + (1:sum(ind{5} > ind{4}(g) & ind{5} < ind{4}(g+1)));
        else
            sweeps{g} = sweeps{g-1}(end) + (1:sum(ind{5} > ind{4}(g)));
        end
    end
end

% Sort data into tree structure
ct = 0; % counter to match experiments to chan* and head* files created by sigTOOL
swIdx = 0; % sweep index (to extract sweep specific information, e.g. sweep start time)
digIdx = []; % Index of digital channels to put them into one display group
for g = 1:numel(groups) % Go through all groups (a new group can be created by the experimenter by clicking "New group" or "New experiment" in the file menu).
    for s = 1:numel(groups{g}) % Go through all series of the current group
        for t = 1:numel(sweeps{series{s}(1)}) % Go through all traces, assuming all sweeps contain the same number of traces  - sweeps are already bundled in sigTOOL chan* variables
            ct = ct + 1;
            % Metainformation contained in header
            data.series(s).traces(t).header = eval(['in.head' num2str(ct)]);
            data.series(s).traces(t).name = data.series(s).traces(t).header.title;
            data.series(s).traces(t).notes = data.series(s).traces(t).header.comment;
            data.series(s).traces(t).unit = data.series(s).traces(t).header.adc.Units;
            data.series(s).traces(t).displayGroup = t;
            data.series(s).traces(t).displayGroupName = data.series(s).traces(t).name;
            if strcmp(data.series(s).traces(t).header.adc.Units,'*1')
                data.series(s).traces(t).displayGroupSep = 1.2;
                digIdx(t) = 1;
            else
                data.series(s).traces(t).displayGroupSep = 0.5*abs(diff(data.series(s).traces(t).header.adc.YLim));
                digIdx(t) = 0;
            end
            f = 1; % Factor for unit conversion
            % Convert all voltages to mV
            if strfind(data.series(s).traces(t).header.adc.Units, 'V')
                data.series(s).traces(t).unit = 'mV';
                switch data.series(s).traces(t).header.adc.Units
                    case 'V'
                        f = 1000;
                    case 'mV'
                        f = 1;
                    case 'µV'
                        f = 0.001;
                end
            end
            % Convert all currents to nA
            if strfind(data.series(s).traces(t).header.adc.Units, 'A')
                data.series(s).traces(t).unit = 'nA';
                switch data.series(s).traces(t).header.adc.Units
                    case 'A'
                        f = 1e9;
                    case 'mA'
                        f = 1e6;
                    case 'µA'
                        f = 1e3;
                    case 'nA'
                        f = 1;
                    case 'pA'
                        f = 0.001;
                end
            end
            % Actual data, rescaled from integer to actual floating
            % digit value
            data.series(s).traces(t).values = f*(double(eval(['in.chan' num2str(ct) '.adc']))'*data.series(s).traces(t).header.adc.Scale+data.series(s).traces(t).header.adc.DC);
            for sw = 1:size(data.series(s).traces(t).values,1)
                data.series(s).traces(t).trialNames{sw} = ['Sweep #' num2str(sw,'%03d')];
                data.series(s).traces(t).trialNotes{sw} = '';
                swIdx = swIdx + 1;
                data.series(s).traces(t).absoluteStartTime{sw} =  in.FileSource.header{ind{4}(swIdx),4}.SwTimeMATLAB;
            end
            % Reset sweep index if more traces from the same series are
            % left in the loop, so each trace gets the same time
            % information.
            if t < numel(sweeps{series{s}(1)})
                swIdx = swIdx - size(data.series(s).traces(t).values,1);
            end
            % Time information of trace (i.e. channel)
            switch in.FileSource.header{ind{5}(t),5}.TrXUnit  % Make sure that the time vector will be generated with seconds as the unit
                case 's'
                    dt = in.FileSource.header{ind{5}(ct),5}.TrXInterval;
                case 'ms'
                    dt = in.FileSource.header{ind{5}(ct),5}.TrXInterval/1e3;
                case 'µs'
                    dt = in.FileSource.header{ind{5}(ct),5}.TrXInterval/1e6;
            end
            data.series(s).traces(t).timeVector = [in.FileSource.header{ind{5}(ct),5}.TrXStart dt];
        end
        % Comments and labels of series
        data.series(s).expNotes = in.FileSource.header{ind{3}(s),3}.SeComment; % Not contained in header
        data.series(s).expName = in.FileSource.header{ind{3}(s),3}.SeLabel; % Also contained in each header (head*.Group.Label)
        data.series(s).name = [pn filesep fName];
        anCt = 0;
        if ~isempty(digIdx)
            for t = 1:numel(sweeps{series{s}(1)})
                if digIdx(t)
                    data.series(s).traces(t).displayGroup = sum(~digIdx)+1;
                    data.series(s).traces(t).displayGroupName = 'Dig In';
                    data.series(s).traces(t).values = scaleMinMax(data.series(s).traces(t).values);
                else
                    anCt = anCt + 1;
                    data.series(s).traces(t).displayGroup = anCt;
                end
            end
        end
    end
    
end
data.sourceFile = in.head1.source;

% Save for later use
save([fn '.mat'], 'data','groups','series','sweeps','-v7.3');
if nargout
    varargout{1} = [fn '.mat'];
end

