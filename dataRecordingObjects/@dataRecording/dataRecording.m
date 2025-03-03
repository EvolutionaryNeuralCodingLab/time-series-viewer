classdef (Abstract) dataRecording < handle
    properties
        recordingName %(String) The name of the recording
        recordingDir % (String) Full directory containing the recorded session
        dataFileNames % (Cell 1 x N)  array of N recording data file names
        startDate %(1x1) Start date (time) of Recording (matlab date number format)
        endDate %(1x1) End date (time) of Recording (matlab long format)
        samplingFrequency %(1xN) Sampling rate of electrode signals [Hz]
        samplingFrequencyAnalog %(1xN) Sampling rate of analog signals [Hz]
        recordingDuration_ms %(1x1) the total duration of the recording in [ms]
        channelNames % (Cell 1xN) a cell array with the N names of the channels
        channelNumbers % (1xN) an array with integer channel numbers (>=1 integer)
        
        triggerNames %the names of trigger channels (not critical)
        analogChannelNumbers % (1xN) the numbers of channels containing non-electrode analog inputs
        analogChannelNames % (1xN) the numbers of channels containing non-electrode analog inputs
        
        chLayoutNumbers %(MxN) The layout of the channel numbers in physical space arranged in an M by N grid
        chLayoutNames %(Cell MxN)The layout of the channel names in physical space arranged in an M by N grid
        electrodePitch % distance between electrodes (not critical)
        chLayoutPositions % (1xN or 2xN or 3xN) array of electrode position in [x or x,y or x,y,z]
        layoutName %the name of the channel layout (electrode type)
        
        convertData2Double = 1; % if data should be converted to double from the original quantization
        ZeroADValue % the digital zero value fir electrode channels
        MicrovoltsPerAD % the digital to analog conversion value for electrode channels
        MicrovoltsPerADAnalog % the digital to analog conversion value for analog channels
        ZeroADValueAnalog% the digital zero value for analog channels
        datatype        % class of data in the recording
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        n2s % a translation between the number of the channel to the serial number of the channel (in the case where all channels are consecutive)
        dspLowCutFrequency % (1x1) Low-pass cutoff frequency in the Neuralynx DSP (in raw data)
        dspHighCutFrequency % (1x1) High-pass cutoff frequency in the Neuralynx DSP (in raw data)
        nRecordings % (1x1) number of recording files
        
        overwriteMetaData = false; %internal parameter - if true recalculates meta data
        metaDataFile
        includeOnlyDigitalDataInTriggers = false;
        channelNumbersOrignal %check if we can remove this - not clear
    end
    
    properties (SetAccess=protected) %these are properties that are not synchronized or loaded from meta files
        multifileMode %(logical 1x1) if multi files were selected %%%for now - false
        folderMode = false; %false forever (yuval)
    end
    
    properties (Constant, Abstract)
        defaultLocalDir %Default directory from which search starts
        signalBits % the quantization of the sampling card
        numberOfCharFromEndToBaseName %the number of characters from the end of the file name to get to the base file name
    end

    methods (Static)
        function concatenateBinaryFormat(binaryFileNames,varargin)
            parseObj = inputParser;
            parseObj.FunctionName='concatenateBinaryFormat';
            addRequired(parseObj,'binaryFileNames',@iscellstr);
            addParameter(parseObj,'targetFileBase',[],@ischar); % in milliseconds
            if numel(varargin)==1
                disp(parseObj.Results);
                return;
            end
            parseObj.parse(binaryFileNames,varargin{:});
            %make parameter structure
            par=parseObj.Results;

            disp('checking that file names are valid');
            if isempty(binaryFileNames)
                fprintf('Files to merge are not provided in the first argument. Aborting...\n');
                return;
            else
                for i=1:numel(binaryFileNames)
                    %Verify that target extension is "bin" or "dat"
                    if ~any(strcmp(binaryFileNames{i}(end-3:end),'.bin'))
                        error('input files should have a ''.bin'' extension');
                    end
                end
            end
            if isempty(par.targetFileBase)
                par.targetFileBase=[cd filesep 'mergedBinary.bin'];
                fprintf('Target file name not entered! Using this: %s\n',par.targetFileBase);
            end

            %create a folder if target folder does not exist
            [folderName,FileName]=fileparts(par.targetFileBase);
            if ~isfolder(folderName)
                fprintf('Notice the folder you entered does not exist!!! Creating new folder: %s\n',folderName);
                mkdir(folderName);
            end

            outMeta=[par.targetFileBase(1:end-4) '_meta.txt'];
            outAnalog=[par.targetFileBase(1:end-4) '_Analog' par.targetFileBase(end-3:end)];
            outTriggers=[par.targetFileBase(1:end-4) '_Triggers.bin'];

            %check if target files exist and if so abort
            if isfile(outMeta) ||...
                    isfile(outAnalog) ||...
                    isfile(par.targetFileBase) ||...
                    isfile(outTriggers)
                error('Target files exist. Please  delete target files and run again!')
            end

            %verify that meta data files are identical
            checkF1 = javaObject('java.io.File', [binaryFileNames{1}(1:end-4) '_meta.txt']);
            for i=2:numel(binaryFileNames)
                checkF2 = javaObject('java.io.File', [binaryFileNames{i}(1:end-4) '_meta.txt']);
                is_equal = javaMethod('contentEquals','org.apache.commons.io.FileUtils',checkF1, checkF2);
                if ~is_equal
                    error('meta data files for the different recordings are not identical. Aborting merge!');
                end
            end

            %Merging triggers.

            %go over triggers and merge time stamps for each
            fidMeta = fopen([binaryFileNames{1}(1:end-4) '_meta.txt']);
            TextAsCells = textscan(fidMeta, '%s = %s','TextType','string','Delimiter',' = ');
            nTriggersCh = str2double(TextAsCells{2}{strcmp(TextAsCells{:,1},"nTriggerChans")});
            nChans = str2double(TextAsCells{2}{strcmp(TextAsCells{:,1},"nChans")});
            fclose(fidMeta);

            %open files for read and write
            fidTargetTrig = fopen(outTriggers, 'w');
            for i=1:numel(binaryFileNames)
                tmpFile=[binaryFileNames{i}(1:end-4) '_Triggers.bin'];
                if exist(tmpFile,'file')
                    fid(i)=fopen(tmpFile,'r');
                    nTriggers(:,i)=double(fread(fid(i),nTriggersCh,'*uint32')-1);
                else
                    error('Some trigger files are missing! Aborting!');
                end

                fidData=fopen(binaryFileNames{i},'r');
                fseek(fidData, 0, 'eof');
                position = ftell(fidData);
                nTotSamplesRecording(i)=floor(position/2/nChans); %2 bytes per sample
                fclose(fidData);
            end
            sampleAddition_unit32=[0 cumsum(nTotSamplesRecording(1:end-1))];

            for i=1:nTriggersCh
                for j=1:numel(binaryFileNames)
                    if nTriggers(i,j)~=0
                        tmp=fread(fid(j),nTriggers(i,j),'*uint32');
                        fwrite(fidTargetTrig,tmp+uint32(sampleAddition_unit32(j)),'*uint32');
                    end
                end
            end

            fclose(fidTargetTrig);
            for i=1:numel(binaryFileNames)
                fclose(fid(i));
            end

            %create new meta data file which is identical to the original files.
            copyfile([binaryFileNames{1}(1:end-4) '_meta.txt'],outMeta);

            %merge analog files and data files
            fprintf('Starting to copy files...\nImportant!!! Matlab command while files are still being copied in the background!\n')
            if isunix
                systemStrData=['cat'];
                systemStrDataAnalog=['cat'];
                for i=1:numel(binaryFileNames)
                    systemStrData=[systemStrData ' ' binaryFileNames{i}];
                    systemStrDataAnalog=[systemStrDataAnalog ' ' binaryFileNames{i}(1:end-4) '_Analog' binaryFileNames{i}(end-3:end)];
                end
                [status,cmdout] = system([systemStrData ' > ' par.targetFileBase],'-echo');
                [status,cmdout] = system([systemStrDataAnalog ' > ' outAnalog],'-echo');

            elseif IsWindows
                systemStrData=['copy ' binaryFileNames{1} '/b'];
                systemStrDataAnalog=['copy ' binaryFileNames{1}(1:end-4) '_Analog' binaryFileNames{1}(end-3:end) '/b'];
                for i=2:numel(binaryFileNames)
                    systemStrData=[systemStrData '+' binaryFileNames{i} '/b'];
                    systemStrDataAnalog=[systemStrDataAnalog '+' binaryFileNames{i}(1:end-4) '_Analog' binaryFileNames{i}(end-3:end) '/b'];
                end
                [status,cmdout] = system([systemStrData ' ' par.targetFileBase],'-echo');
                [status,cmdout] = system([systemStrDataAnalog ' ' outAnalog],'-echo');
            end

        end
    end

    methods
        function delete(obj) %closing all open files when object is deleted
            obj=closeOpenFiles(obj);
        end
        function obj=closeOpenFiles(obj)
        end
        function [V_uV,t_ms]=getData(obj,channels,startTime_ms,window_ms,name)
            %Extract recording data from file to memory
            %Usage: [V_uV,t_ms]=obj.getData(channels,startTime_ms,window_ms);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the recording (if empty takes the default name)
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        t_ms - A time vector relative to recording start (t=0 at start)
        end
        function [V_uV,T_ms]=getAnalogData(obj,channels,startTime_ms,window_ms,name)
            %Extract recording data from file to memory
            %Usage: [V_uV,t_ms]=obj.getAnalogData(channels,startTime_ms,window_ms);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the recording (if empty takes the default name)
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        t_ms - A time vector relative to recording start (t=0 at start)
        end
        function [T_ms]=getTrigger(obj,startTime_ms,window_ms,name)
            %Extract triggers from file Neuralynx recording
            %Usage : [T_ms]=obj.getTrigger(startTime_ms,endTime_ms,TTLbits)
            %Input : startTime_ms - start time [ms].
            %        window_ms - the window duration [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        name - the name of the trigger (if empty takes the default name)
            %Output: T_ms - trigger times [ms]
        end
        function [D,T_ms]=getDigitalData(obj,startTime_ms,window_ms,name)
            %Extract MCRack digital data from file to memory
            %Usage: [D,T_ms]=getDigitalData(startTime_ms,window_ms,name)
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the stream (if not entered, default name is used)
            %Output: D - A 3D matrix [nChannels x nTrials x nSamples] with digitalData waveforms across specified channels and trials
            %        T_ms - A time vector relative to recording start (t=0 at start)
        end
        function saveMetaData(obj)
            %Save object properties (metaData) to file
            %Usage : obj.saveMetaData;
            %Input : []
            props.metaClassData=metaclass(obj);
            props.allPropName={props.metaClassData.PropertyList.Name}';
            props.allPropIsConstant=cell2mat({props.metaClassData.PropertyList.Constant}');
            props.allPropSetAccess={props.metaClassData.PropertyList.SetAccess}';
            
            pNonConstantProps=find(~props.allPropIsConstant & ~strcmp(props.allPropSetAccess,'protected'));
            for i=1:numel(pNonConstantProps)
                metaData.(props.allPropName{pNonConstantProps(i)})=obj.(props.allPropName{pNonConstantProps(i)});
            end
            if iscell(obj.recordingDir)
                recordingDir=obj.recordingDir{1};
                recordingName = strjoin(cellfun(@(x) x(1:end-3), obj.dataFileNames,'UniformOutput',false),'-');
            else
                recordingDir=obj.recordingDir;
                recordingName = obj.recordingName; %make sure this always work
            end
            save(obj.metaDataFile,'metaData');
        end
        
        function [X,Y,Z]=getElectrodePositions(obj,electrodePitch)
            %if recording object contains electrode positions, use these, if not
            if ~isempty(obj.chLayoutPositions)
                disp('Getting positions from layout files');
                X=obj.chLayoutPositions(1,:);
                Y=obj.chLayoutPositions(2,:);
            else
                if nargin==2
                    obj.electrodePitch=electrodePitch;
                elseif nargin==1 & isempty(obj.electrodePitch)
                    obj.electrodePitch=100;
                end
                disp(['Getting positions from grid layout, assuming pitch of ' num2str(obj.electrodePitch) 'um !!!!!']);
                
                %Build inverse map between electrode and location
                [meshX,meshY]=meshgrid(1:size(obj.chLayoutNumbers,1),1:size(obj.chLayoutNumbers,2));
                X(obj.chLayoutNumbers(~isnan(obj.chLayoutNumbers)))=meshX(~isnan(obj.chLayoutNumbers))*obj.electrodePitch;
                Y(obj.chLayoutNumbers(~isnan(obj.chLayoutNumbers)))=meshY(~isnan(obj.chLayoutNumbers))*obj.electrodePitch;
            end
            Z=zeros(size(Y));
        end
        
        function deleteMetaData(obj)
            if ~iscell(obj.recordingDir)
                delete(obj.metaDataFile);
            else
                for i=1:numel(obj.recordingDir)
                    delete(obj.metaDataFile);
                end
            end
        end
        
        function obj=loadMetaData(obj,fileName)
            %Load object properties (metaData) from file
            %Usage : obj.loadMetaData;
            %Input : fileName - if entered loads meta data from this file, else loads data from main data directory
            try
                oldRecordingDir=obj.recordingDir;
                if ~iscell(obj.recordingDir) %regular recording
                    if nargin==2
                        load(fileName);
                    else
                        load(obj.metaDataFile,'metaData');
                    end
                    fieldNames=fieldnames(metaData);
                    for i=1:numel(fieldNames)
                        obj.(fieldNames{i})=metaData.(fieldNames{i});
                    end
                else %multi file recording
                    recordingName = strjoin(cellfun(@(x) x(1:end-3), obj.dataFileNames,'UniformOutput',false),'-');
                    for i=1:numel(obj.recordingDir)
                        if nargin==2
                            load(fileName{i});
                        else
                            load(obj.metaDataFile,'metaData');
                        end
                        fieldNames=fieldnames(metaData);
                        for j=1:numel(fieldNames)
                            if numel(metaData.(fieldNames{j}))==1
                                obj.(fieldNames{j})(i)=metaData.(fieldNames{j});
                            else
                                if ~iscell(obj.(fieldNames{j}))==1
                                    obj.(fieldNames{j})=cell(1,numel(obj.recordingDir));
                                end
                                obj.(fieldNames{j}){i}=metaData.(fieldNames{j});
                            end
                        end
                    end
                end
                obj.recordingDir=oldRecordingDir;
            catch errorMsg
                disp('Error while extracting fields from meta data. Trying re-extract meta data...');
                obj=extractMetaData(obj);
            end
        end
        
        function obj=loadChLayout(obj,layoutName)
            if nargin==1
                %checks for a .chMap file with the recording name in the same folder of the recording and extract the layout information
                %txt should correspond to layout file name on path
                if iscell(obj.recordingDir)
                    recordingDir=obj.recordingDir{1};
                else recordingDir=obj.recordingDir;
                end
                
                chMapFiles=dir([recordingDir filesep '*.chMap']);
                chMapFiles={chMapFiles.name};
                if numel(chMapFiles)>1
                    disp('Found more than one channel map files!!!!! Using one of them');
                    chMapFiles=chMapFiles(1);
                end
                
                switch numel(chMapFiles)
                    case 0 %not channel map file found
                        disp('No .chMap files were found for this recording');
                        return;
                    case 1 %there is only one channel map file, this file will apply to all the recordings
                        chMapFiles=[recordingDir filesep chMapFiles{1}];
                    otherwise %there are several files, in which case each recording should have its own channel map file with the appropriate name
                        chMapFiles=dir([recordingDir filesep obj.recordingName(1:end-numel(obj.fileExtension)-1) '*.chMap']);
                        chMapFiles={chMapFiles.name};
                        if numel(chMapFiles)~=1
                            disp('Channel map file name (*.chMap) does not correspond to the recording file name');
                            return;
                        else
                            chMapFiles=[recordingDir filesep chMapFiles{1}];
                        end
                end
                
                layoutName = importdata(chMapFiles);
                if isempty(layoutName)
                    error('Channel layout was not extracted successfully from channel map file (*.chMap) ,check file name and content!');
                end
                fprintf('Using the file %s for extracting channel map\n',chMapFiles);
            else
                if ~isempty(layoutName)
                %Convert to cell array
                if ~iscell(layoutName)
                    layoutName={layoutName};
                end
                %remove the layout prefix if entered
                if all(layoutName{1}(1:7)=='layout_')
                    layoutName{1}=layoutName{1}(8:end);
                end
                else
                    warning('Provided layout name is empty! Layout not updated!')
                    return;
                end
            end
            
            
            try
                obj.layoutName=layoutName{1};
                allElectrodes=regexp(layoutName{1},',','split');lastElectrode=0;
                if numel(allElectrodes)>1
                    disp('Multiple electrode detected, stacking them horizontally in display');
                end
                
                obj.chLayoutNumbers=[];
                obj.chLayoutNames={};
                obj.chLayoutPositions=[];
                lastElectrode=0;
                for i=1:numel(allElectrodes)
                    elecString=regexp(allElectrodes{i},'_','split');
                    electrodePitch(i)=str2num(elecString{1});
                    
                    load(['layout_' allElectrodes{i} '.mat']);
                    chLayoutNumbers=[obj.chLayoutNumbers;En+lastElectrode];
                    chLayoutNames=[obj.chLayoutNames;Ena];
                    chLayoutPositions=[obj.chLayoutPositions,[Enp(1,:);Enp(2,:)]]; %check why this was changed.
                    %chLayoutPositions=[obj.chLayoutPositions,[Enp(1,:);Enp(2,:)+max(Enp(2,:))+200]];

                    lastElectrode=size(obj.chLayoutPositions,2);
                    fprintf('Channel map with pitch %d and layout %s extracted\n',obj.electrodePitch,elecString{2});
                end
                %update class at the end in case there was an error and a previous layout already existed
                obj.chLayoutNames=chLayoutNames;
                obj.electrodePitch=electrodePitch;
                obj.chLayoutNumbers=chLayoutNumbers;
                obj.chLayoutPositions=chLayoutPositions;
                
                %check that all recorded channels are contained within the layout
                if numel(obj.channelNumbers)>numel(intersect(obj.channelNumbers,En(:)))
                    warning(['Notice that some of the recorded channels (%d channels) are not contained in the layout file (%d channels), this may result in errors in further analysis!...' ...
                        '\nIf more than one probe is recorded, probe names can be comma seperate in chMap file'],numel(obj.channelNumbers),sum(~isnan(En(:))));
                end
            catch
                fprintf('Failed to extract channel map!!!! check that the name was entered correctly\n');
            end
            
        end
        
        function []=getKiloSort(obj,tempFilesFolder,varargin)
            if nargin<2
                fprintf('First argument should be a temporary folder. It''s better if its local.\nThis folder is important for later sorting with phy.\n');
            end
            parseObj = inputParser;
            parseObj.FunctionName='getKiloSort';
            addRequired(parseObj,'tempFilesFolder',@ischar);
            addParameter(parseObj,'tStart',0,@isnumeric); % in milliseconds
            addParameter(parseObj,'tEnd',Inf,@isnumeric); % in milliseconds
            addParameter(parseObj,'correctDrift',1,@isnumeric); % if to run drift correction (in any case drift will be examined)
            addParameter(parseObj,'ops',[]); % a structure with parameters
            addParameter(parseObj,'runManualCurationPhy',1); %run manual curation in phy after sorting
            addParameter(parseObj,'overwrite',0,@isnumeric); %overwrite everything
            addParameter(parseObj,'runOnlySpikeDetection',0,@isnumeric); %stops after spike detection. Can be run again later
            addParameter(parseObj,'saveOnlyPreClusteringData',0,@isnumeric); %save the data just before clustering so that this step could be repeated
            addParameter(parseObj,'loadPreClustering',0,@isnumeric); % load pre clustering data.

            addParameter(parseObj,'outFolder',fullfile(obj.recordingDir,['kiloSortResults_',obj.recordingName]),@ischar);
            if numel(varargin)==1
                disp(parseObj.Results);
                return;
            end
            parseObj.parse(tempFilesFolder,varargin{:});
            %make parameter structure
            par=parseObj.Results;
            
            % generate a recording specific identifier
            recID=num2str(sum(par.outFolder));
            
            % generate config structure
            rootH = par.tempFilesFolder; %where to save temp files for spike sorting (should be a fast drive)
            
            %updsampingFactor=1.5;
            %(pre-processing)
            ops.trange    = [par.tStart',par.tEnd']/1000; % (pre-processing) time range to sort - move to time units of seconds
            ops.NchanTOT  = numel(obj.channelNumbers); % (pre-processing) total number of channels in your recording
            ops.fproc   = fullfile(rootH,['temp_wh_' recID '.dat']); % proc file on a fast SSD
            ops.fbinary = fullfile(obj.recordingDir, obj.dataFileNames{1});

            ops.fs = obj.samplingFrequency(1)*1.5;% sample rate
            ops.fshigh = 300;% frequency for high pass filtering (150)
            ops.fslow = 3000; % frequency for low pass filtering (optional)
            ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry
            ops.parfor              = 1; % whether to use parfor to accelerate some parts of the algorithm

            ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
            ops.NT                  = 64*1024+ ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory).
            ops.whiteningRange      = Inf; % number of channels to use for whitening each channel (defualt 32)
            ops.nSkipCov            = 25; % compute whitening matrix from every N-th batch
            ops.scaleproc           = 200;   % int16 scaling of whitened data
            ops.useRAM              = 0; % not yet available

            %(extract_spikes)
            ops.spkTh           = -6;      % spike threshold in standard deviations (-6)
            ops.long_range = [30 6]; % ranges to detect isolated peaks ([30 6])

            %(extract_spikes)
            ops.Th = [25 20];% Thresholds on spike detection used during the optimization Th(1) or during the final pass Th(2). These thresholds are applied to the template projections, not to the voltage. Typically, Th(1) is high enough that the algorithm only picks up sortable units, while Th(2) is low enough that it can pick all of the spikes of these units. It doesn't matter if the final pass also collects noise: an additional per neuron threshold is set afterwards, and a splitting step ensures clusters with multiple units get split.

            %ops.Nfilt               =
            % 1024; % max number of clusters
            ops.nNeigh              = 32; %For visualization only - number of neighboring templates to retain projections (default 16).
            ops.nfilt_factor        = 32; % max number of clusters per good channel (even temporary ones)
            ops.minfr_goodchannels = 0; % minimum firing rate on a "good" channel (0 to skip)
            ops.lam = 10;% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot).  The individual spike amplitudes are biased towards the mean of the cluster by this factor; 50 is a lot, 0 is no bias.
            % large lam means amplitudes are forced around the mean ([10 30 30]) // The higher these values are, the more Kilosort ensures spikes within a cluster have close amplitudes
            ops.AUCsplit = 0.3;% splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1). Threshold on the area under the curve (AUC) criterion for performing a split in the final step. If the AUC of the split is higher than this, that split is considered good. However, a good split only goes through if, additionally, the cross-correlogram of the split units does not contain a big dip at time 0.
            ops.minFR = 1/500;% minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
            ops.momentum = [20 400];% number of samples to average over (annealed from first to second value)
            ops.sigmaMask = 130;% spatial constant in um for computing residual variance of spike (default 30)
            ops.ThPre = 8;% threshold crossings for pre-clustering (in PCA projection space) (originally 8)
            ops.CAR=1; %median filter over channels (seems to slightly improve the spike quality).
            ops.sig= 10;  % spatial smoothness constant for registration (20)
            ops.nblocks    = 5; % blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option.
            %ops.criterionNoiseChannels = 0.1; % fraction of "noise" templates allowed to span all channel groups (see createChannelMapFile for more info).

            % ops.reorder         = 1;       % whether to reorder batches for drift correction. (changing the settings below can lead to fatal errors!!!)
            ops.nskip           = 25;  % how many batches to skip for determining spike PCs (changing the settings below can lead to fatal errors!!!)
            ops.nPCs                = 3; % how many PCs to project the spikes into (changing the settings below can lead to fatal errors!!!)
            if ~isempty(par.ops)
                fields = fieldnames(par.ops);
                for i = 1:length(fields)
                    ops.(fields{i}) = par.ops.(fields{i});
                end
            end

            if ~strcmp(class(obj),'binaryRecording') %check if this is a binary recording.
                fprintf('\nKilosort can only run on binary files, use the export2Binary method to first convert the data.\n Then switch to the binaryRecording object and run again\n');return;
            end
            
            fprintf('\nConverting layout to kilosort format');
            obj.convertLayoutKSort;
            % is there a channel map file in this folder?
            channelMapFile = dir(fullfile(obj.recordingDir, 'chan*.mat'));
            if ~isempty(obj.recordingDir)
                ops.chanMap = fullfile(obj.recordingDir, channelMapFile.name);
            else
                fprintf('\nThe Channel map is empty or missing. Please provide and run again!');return;
            end
            
            [kilosortPath]=which('kilosort');
            if isempty(kilosortPath)
                fprintf('Kilosort was not found, please add it to the matlab path and run again (e.g. addpath(genpath(''/home/mark/Documents/MATLAB/Kilosort''))');
                return;
                %addpath(genpath('/home/mark/Documents/MATLAB/Kilosort')) % path to kilosort folder
            end

            writeNPYPath=which('writeNPY.m');
            rez2PhyPath2=which('rezToPhy2.m');
            if isempty(rez2PhyPath2) || isempty(writeNPYPath)
                fprintf('\nrez2Phy2 (kilosort3) or writeNPY was not found, please add it to the matlab path and run again\n');
                return;
            end
            %ch2Remove=[18 22 23 30 31]

            [~,recFolder]=fileparts(obj.recordingDir);
            expName=['kilosortRez_' recFolder];
            tmpSaveFile=[rootH filesep expName '_' num2str(sum(par.outFolder)) '.mat']; %create a unique name for every experiment

            %loads previously processes data
            if ~isfolder(rootH)
                mkdir(rootH);
            end
            if isfile(tmpSaveFile) & ~par.overwrite
                load(tmpSaveFile);
            end

            if ~par.loadPreClustering
                %check existance of preprocessing in previously saved data
                if ~exist('rezPreProc','var') || par.overwrite==true
                    rezPreProc = preprocessDataSub(ops);
                    save(tmpSaveFile,'rezPreProc','-v7.3');
                else
                    disp('pre processing loaded from last saved version');
                end

                if ~exist('rezShift','var') || par.overwrite==true
                    rezShift                = datashift2(rezPreProc,par.correctDrift);
                    save(tmpSaveFile,'rezShift','-append');
                else
                    disp('Shift processing loaded from last saved version');
                end

                %use kilosort 3
                if ~exist('rezSpk','var') || par.overwrite==true
                    %Adding a try catch inside extract_spikes arround st(5,:) = cF; in
                    %the function extract_spikes in the folder clustering of
                    %kilosort seems to solve the problem of crashing during execution
                    [rezSpk, st3, tF]     = extract_spikes(rezShift);
                    save(tmpSaveFile,'rezSpk','st3','tF','-append');
                else
                    disp('Spike extraction loaded from last saved version');
                end

                if ~isfolder(par.outFolder)
                    mkdir(par.outFolder);
                end

                if par.runOnlySpikeDetection
                    return;
                end

                if par.saveOnlyPreClusteringData
                    disp('A copy of the pre-clustering results will not be deleted after sorting. Make sure to delete manually when not needed!')
                    save([par.outFolder,filesep,'preClusteringResults.mat'],'rezSpk','tF','st3'); %the parameter file ops is included
                end
            else
                disp(['Loading pre clustered results from: ' par.outFolder,filesep,'preClusteringResults.mat'])
                load([par.outFolder,filesep,'preClusteringResults.mat']); %the parameter file ops is included
                %update the field names in the ops file saved in rezSpk
                for i=fieldnames(ops)'
                    rezSpk.ops.(i{1})=ops.(i{1});
                end
            end

            if ~exist('rez','var') || par.overwrite || par.loadPreClustering
                try
                    rez                = template_learning(rezSpk, tF, st3);
                catch
                    disp('Error with CUDA, trying again!');
                    rez                = template_learning2(rezSpk, tF, st3);
                    %pause(2);
                    %rez                = template_learning(rezSpk, tF, st3);
                end
                [rez, st3, tF]     = trackAndSort(rez);
                rez                = final_clustering(rez, tF, st3);
                rez                = find_merges(rez, 1);

                % correct times for the deleted batches
                %rez=correct_time(rez);
                if ~par.loadPreClustering
                    save(tmpSaveFile,'rez','-append');
                end
            else
                disp('Final template learning and clustering loaded from last saved version');
            end

            fprintf('Done kilosort\nSaving results and exporting Phy templates to %s\n',par.outFolder);
            if isfolder(par.outFolder)
                fprintf('Deleting previous kilosort results...\n');
                delete([par.outFolder filesep '*.npy']);
                delete([par.outFolder filesep '*.tsv']);
                delete([par.outFolder filesep 'params.py']);
                delete([par.outFolder filesep 'phy.log']);
            end
            rezToPhy2(rez, par.outFolder);

            %save settings file
            disp(['Saving options file in ' par.outFolder filesep 'ops.mat']);
            save([par.outFolder filesep 'ops.mat'],'ops');
            
            if ~par.loadPreClustering
                delete(tmpSaveFile);
            end

            %run manual curation with Phy
            if par.runManualCurationPhy
                %conda('activate','phy2')
                setenv('PATH','/home/mark/anaconda3/envs/phy2/bin:/sbin:/bin:/usr/bin:/usr/local/bin:/snap/bin:/home/mark/anaconda3/bin');%set phy2 environment
                system(['phy template-gui "' fullfile(obj.recordingDir,['kiloSortResults_',obj.recordingName]) filesep 'params.py"']);%run phy on current analysis object
            end
        end

        function [qMetric,unitType] = getBombCell(obj,pathToKSresults,GUIbc,rerun)%GUIbc,rerun)

            if nargin <3 
                GUIbc = 0;
            end

            if nargin <4 
                rerun = 0;
            end

            %%%%%Set paths

            %These paths below are the paths you will need to input to load data and 
            % save the computed quality metrics / ephys properties. 
            % All metrics relating to raw data (amplitude, signal to noise ratio) will not be computed. 

            ephysKilosortPath = pathToKSresults;
%             ephysRawFile = dir([pathToKSresults '*ap*.*bin']); % path to your raw .bin or .dat data
%             ephysMetaDir = dir([pathToKSresults '*ap.meta']); % path to your .meta or .oebin meta file
            savePath = [pathToKSresults 'qMetrics']; % where you want to save the quality metrics

            files = dir(pathToKSresults);
            fileNames = {files.name};
            APbin = char(fileNames(contains(fileNames,'tcat.imec0.ap.bin')));
            METAbin = char(fileNames(contains(fileNames,'tcat.imec0.ap.meta')));
           
            ephysRawFile = dir(fullfile(pathToKSresults,APbin)); % path to yourraw .bin or .dat data
            ephysMetaDir = dir(fullfile(pathToKSresults,METAbin)); % path to your .meta or .oebin meta file

            kilosortVersion = 4; % if using kilosort4, you need to have this value kilosertVersion=4. Otherwise it does not matter.
            gain_to_uV = NaN; % use this if you are not using spikeGLX or openEphys to record your data. this value,
            % when mulitplied by your raw data should convert it to  microvolts.

            %%%%Load data
            %This function loads are your ephys data. Use this function rather than any custom one
            %as it handles zero-indexed values in a particular way. 

            [spikeTimes_samples, spikeClusters, templateWaveforms, templateAmplitudes, pcFeatures, ...
                pcFeatureIdx, channelPositions] = bc.load.loadEphysData(ephysKilosortPath, savePath);

            %%%%Run quality metrics
            %%%%Set your paramaters.
            %These define both how you will run quality metrics and how thresholds will be applied to quality metrics 
            % to classify units into good/MUA/noise/non-axonal. 
            % This function loads default, permissive values. 
            % It's highly recommended for you to iteratively tweak these values to find values that suit your particular use case!

            param = bc.qm.qualityParamValues(ephysMetaDir, ephysRawFile, ephysKilosortPath, gain_to_uV, kilosortVersion);

            %%%%Specific paremeters to tweek for lizards:
            param.maxWvBaselineFraction = 0.5; %default value = 0.3

            param.computeTimeChunks = 1; %default = 0

            param.deltaTimeChunk = 20*60; %default = 360

            param.maxPercSpikesMissing = 30;%default = 20

            param.tauR_valuesMin = 0.002;%default = 0.002

            param.tauR_valuesMax = 0.002;%default = 0.002

            param.hillOrLlobetMethod =1;%default = 1

            param.maxDrift = 200;%default = 100

            param.minPresenceRatio = 0.15;%default = 0.7

            param.removeDuplicateSpikes = 1; %default = 0;

            param.minAmplitude = 20;%default = 20

            param.splitGoodAndMua_NonSomatic = 0; %default = 0; %%%Not working, raise issue

            param.minNumSpikes = 200; %default = 300

            param.maxRPVviolations = 0.3; %fraction, default = 0.1 %%Weirdly, the refractory period violations are not reflected in the autocorrelation. How is it calculated?


            %param.nChannels must correspond to the total number of channels in your raw data, including any sync channels. 
            % For Neuropixels probes, this value should typically be either 384 or 385 channels. 
            % param.nSyncChannels must correspond to the number of sync channels you recorded. This value is typically 1 or 0.

            param.nChannels = 385;
            param.nSyncChannels = 1;

            % if using SpikeGLX, you can use this function:
            if ~isempty(ephysMetaDir)
                if endsWith(ephysMetaDir.name, '.ap.meta') %spikeGLX file-naming convention
                    meta = bc.dependencies.SGLX_readMeta.ReadMeta(ephysMetaDir.name, ephysMetaDir.folder);
                    [AP, ~, SY] = bc.dependencies.SGLX_readMeta.ChannelCountsIM(meta);
                    param.nChannels = AP + SY;
                    param.nSyncChannels = SY;
                end
            end

            %%%%%Run all your quality metrics!
            % This function runs all quality metrics, saves the metrics in your savePath folder and outputs some global summary 
            % plots that can give you a good idea of how things went.

            qMetricsExist = ~isempty(dir(fullfile(savePath, 'qMetric*.mat'))) || ~isempty(dir(fullfile(savePath, 'templates._bc_qMetrics.parquet')));

            if qMetricsExist == 0 || rerun
                [qMetric, unitType] = bc.qm.runAllQualityMetrics(param, spikeTimes_samples, spikeClusters, ...
                templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, savePath);
            else
                [param, qMetric] = bc.load.loadSavedMetrics(savePath);
                unitType = bc.qm.getQualityUnitType(param,qMetric,savePath);
                %unitType = bc_getQualityUnitType(param, sortrows(qMetric,'maxChannels'), savePath);
            end

            %%%%%Inspect
            % After running quality metrics, espacially the first few times, it's a good idea to inspect your data and the quality metrics using the built-in GUI. Use your keyboard to navigate the GUI:
            % left/right arrow : toggle between units
            % u  : brings up a input dialog to enter the unit you want to go to
            % g  : go to next good unit
            % m : go to next multi-unit
            % n  : go to next noise unit
            % a  : go to next non-somatic unit ("a" is for axonal)
            % up/down arrow : toggle between time chunks in the raw data

            if GUIbc == 1
                bc.load.loadMetricsForGUI;

                loadRawTraces = 0;
                unitQualityGuiHandle = bc.viz.unitQualityGUI_synced(memMapData, ephysData, qMetric, forGUI, rawWaveforms, ...
                    param, probeLocation, unitType, loadRawTraces);
            end
           
% Previous BC run code 
%             %%%%%% 
% 
%             files = dir(pathToKSresults);
% 
%             fileNames = {files.name};
% 
%             APbin = char(fileNames(contains(fileNames,'tcat.imec0.ap.bin')));
%             METAbin = char(fileNames(contains(fileNames,'tcat.imec0.ap.meta')));
%             
%             ephysRawDir = dir(fullfile(pathToKSresults,APbin)); % path to yourraw .bin or .dat data
%             ephysMetaDir = dir(fullfile(pathToKSresults,METAbin)); % path to your .meta or .oebin meta file
%             saveLocation = pathToKSresults(1:strfind(pathToKSresults,'catgt')-2);
%             savePath = fullfile(saveLocation, 'qMetrics');
%             decompressDataLocal = saveLocation; % where to save raw decompressed ephys data
% 
% 
%             %%% load data
%             [spikeTimes_samples, spikeTemplates, templateWaveforms, templateAmplitudes, pcFeatures, ...
%                 pcFeatureIdx, channelPositions] = bc_loadEphysData(pathToKSresults);
% 
%             %%% detect whether data is compressed, decompress locally if necessary
%             rawFile = bc_manageDataCompression(ephysRawDir, decompressDataLocal);
% 
%             %%% which quality metric parameters to extract and thresholds
%             param = bc_qualityParamValues(ephysMetaDir, rawFile, pathToKSresults); %for unitmatch, run this:
%             % param = bc_qualityParamValuesForUnitMatch(ephysMetaDir, rawFile, ephysKilosortPath, gain_to_uV)
% 
% %             param.firstPeakRatio = 1.3;
% %             param.minThreshDetectPeaksTroughs = 0.1;
% %             param.minThreshDetectPeaksTroughs = 
% % 
% %             param.minSpatialDecaySlope
% 
% 
%            % param.maxWvBaselineFraction = 
% 
%             %%% compute quality metrics
%           
%             qMetricsExist = ~isempty(dir(fullfile(savePath, 'qMetric*.mat'))) || ~isempty(dir(fullfile(savePath, 'templates._bc_qMetrics.parquet')));
% 
%             if qMetricsExist == 0 || rerun
%                 [qMetric, unitType] = bc_runAllQualityMetrics(param, spikeTimes_samples, spikeTemplates, ...
%                     templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, savePath);
%             else
%                 [param, qMetric] = bc_loadSavedMetrics(savePath);
%                 unitType = bc_getQualityUnitType(param, sortrows(qMetric,'maxChannels'), savePath);
%             end
% 
% 
%             %%% view units + quality metrics in GUI
%             % load data for GUI
% 
%             if GUIbc == 1
%                 loadRawTraces = 0; % default: don't load in raw data (this makes the GUI significantly faster)
%                 bc_loadMetricsForGUI;
% 
%                 % GUI guide:
%                 % left/right arrow: toggle between units
%                 % g : go to next good unit
%                 % m : go to next multi-unit
%                 % n : go to next noise unit
%                 % up/down arrow: toggle between time chunks in the raw data
%                 % u: brings up a input dialog to enter the unit you want to go to
%                 unitQualityGuiHandle = bc_unitQualityGUI(memMapData, ephysData, qMetric, forGUI, rawWaveforms, ...
%                     param, probeLocation, unitType, loadRawTraces);
%             end
% 
%             %%% Save figures
% 
%             if qMetricsExist == 0 || rerun
%                 cd(savePath)
%                 savefig(figure(1), 'templateWaveform.fig')
%                 savefig(figure(2), 'qualityMetrics.fig')
%             end

        end

        function [spkData]=convertPhySorting2tIc(obj,pathToPhyResults,tStart,BombCelled)
            spkData=[];
            if nargin==1
                pathToPhyResults=fullfile(obj.recordingDir,['kiloSortResults_',obj.recordingName]);
                fprintf('Sorting results path not provided, using this path:\n%s\n',pathToPhyResults);
            else
                if isempty(pathToPhyResults)
                    pathToPhyResults=fullfile(obj.recordingDir,['kiloSortResults_',obj.recordingName]);
                    fprintf('Sorting results path not provided, using this path:\n%s\n',pathToPhyResults);
                end
            end
            if ~isfile([pathToPhyResults filesep 'cluster_info.tsv'])
                fprintf('Manual annotation phase with phy not completed! Please run again after completing!!!\n');
            end
            if ~isfile([pathToPhyResults filesep 'cluster_info.tsv'])
                disp('cluster_info.tsv not found! Please first perform manual annotation using phy and try again.');
                return;
            end
            
            saveFileAll=[pathToPhyResults filesep 'sorting_tIc_All.mat'];
            saveFileValid=[pathToPhyResults filesep 'sorting_tIc.mat'];
            
            if nargout==1
                if isfile(saveFileValid)
                    spkData=load(saveFileValid);
                    return;
                end
            end

            if nargin <3
                BombCelled = 0;
            end
            
            readNPYPath=which('readNPY.m');
            if isempty(readNPYPath)
                fprintf('readNPY was not found, trying to add to path please add it to the matlab path and run again\n');
                return;
            end
            %{
            general info:
            template_id - ranges from 0 n_templates-1
            each template is the spike signature on relevant electrodes - spike template multiplied by amplitude for every spike should give ~spike shape
            cluster - a set of spike supposedly fired by the same neuron
            cluster have a cluster_id - 0 - n_clusters-1
            when spikes are removed or added to a cluster the cluster id changes
            This means that we can not use clusters to extract spike templates. We need to regenerate them
            %}
            
            clusterTable=readtable([pathToPhyResults filesep 'cluster_info.tsv'],'FileType','delimitedtext');
            clusterTable=sortrows(clusterTable,'ch');
            spike_clusters = readNPY([pathToPhyResults filesep 'spike_clusters.npy']);
            %spike_templates = readNPY([pathToPhyResults filesep 'spike_templates.npy']);
            spike_times = readNPY([pathToPhyResults filesep 'spike_times.npy']);
            %labelKS = clusterTable.KSLabel;
            clusterGroup=readtable([pathToPhyResults filesep 'cluster_group.tsv'],'FileType','delimitedtext');
            
            if BombCelled
                label = clusterTable.bc_unitType;
            else
                label = clusterTable.group;
            end
            nSpks=clusterTable.n_spikes;
            neuronAmp=clusterTable.amp;
            %spikeShapes=readNPY([pathToPhyResults filesep 'templates.npy']); %check if this needs to be sorted
            %check for clusters on the same electrode
            [uab,a,b]=unique(clusterTable.ch);
            ic=zeros(4,numel(uab));
            t=cell(1,numel(uab));
            currentIdx=0;prevCh=-1;
            nClusters=numel(clusterTable.ch);
            for i=1:nClusters
                if any("id" == string(clusterTable.Properties.VariableNames))
                    t{i}=spike_times(spike_clusters==clusterTable.id(i))'; %changed from id to cluster_id
                else
                    t{i}=spike_times(spike_clusters==clusterTable.cluster_id(i))'; %changed from id to cluster_id
                end
                ic(1,i)=clusterTable.ch(i);
                ic(3,i)=currentIdx+1;
                ic(4,i)=currentIdx+numel(t{i});
                if prevCh==ic(1,i)
                    ic(2,i)=ic(2,i-1)+1;
                else
                    ic(2,i)=1;
                end
                prevCh=ic(1,i);
                currentIdx=ic(4,i);
            end
            t=double(cell2mat(t))/(obj.samplingFrequency(1)/1000);

            if nargin==3
                t=t+tStart;
            end
            if min(ic(1,:))==0
                ic(1,:)=ic(1,:)+1;
            end
            fprintf('Saving results to %s\n',saveFileValid);
            save(saveFileAll,'t','ic','label','neuronAmp','nSpks'); %save full spikes including noise

            if BombCelled
                pValid=strcmp(label,'GOOD')|strcmp(label,'MUA')|strcmp(label,'NON-SOMA');

            else
                pValid=strcmp(label,'good')|strcmp(label,'mua');

            end
            label=label(pValid);
            neuronAmp=neuronAmp(pValid);
            nSpks=nSpks(pValid);
            [t,ic]=RemainNeurons(t,ic,ic(1:2,pValid));
            save(saveFileValid,'t','ic','label','neuronAmp','nSpks');
            
            if nargout==1 %if output is needed and calculation was needed (no saved file existing).
                spkData=load(saveFileValid);
            end
        end

        function [spkData]=convertBCSorting2tIc(obj,pathToBCResults,tStart)
            spkData=[];
            if nargin==1
                pathToBCResults=fullfile(obj.recordingDir,['kiloSortResults_',obj.recordingName]);
                fprintf('Sorting results path not provided, using this path:\n%s\n',pathToBCResults);
            else
                if isempty(pathToBCResults)
                    pathToBCResults=fullfile(obj.recordingDir,['kiloSortResults_',obj.recordingName]);
                    fprintf('Sorting results path not provided, using this path:\n%s\n',pathToBCResults);
                end
            end
            
            saveFileAll=[pathToBCResults filesep 'BCsorting_tIc_All.mat'];
            saveFileValid=[pathToBCResults filesep 'BCsorting_tIc.mat'];
            
            if nargout==1
                if isfile(saveFileValid)
                    spkData=load(saveFileValid);
                    return;
                end
            end
            
            readNPYPath=which('readNPY.m');
            if isempty(readNPYPath)
                fprintf('readNPY was not found, trying to add to path please add it to the matlab path and run again\n');
                return;
            end
            %{
            general info:
            template_id - ranges from 0 n_templates-1
            each template is the spike signature on relevant electrodes - spike template multiplied by amplitude for every spike should give ~spike shape
            cluster - a set of spike supposedly fired by the same neuron
            cluster have a cluster_id - 0 - n_clusters-1
            when spikes are removed or added to a cluster the cluster id changes
            This means that we can not use clusters to extract spike templates. We need to regenerate them
            %}
            
            [spike_times, spike_clusters, templateWaveforms, templateAmplitudes, pcFeatures, ...
                pcFeatureIdx, channelPositions] = bc_loadEphysData(pathToBCResults);


            savePath = pathToBCResults(1:strfind(pathToBCResults,'catgt')-2);
            savePath = fullfile(savePath, 'qMetrics');
            cd(savePath);

            [param, qMetric] = bc_loadSavedMetrics(savePath);
            unitType = bc_getQualityUnitType(param, sortrows(qMetric,'maxChannels'), savePath);
            dspikes = readNPY('spikes._bc_duplicateSpikes.npy'); %Spike is duplicate?

            cd(pathToBCResults);

            %clusterTable=readtable([pathToBCResults filesep 'cluster_info.tsv'],'FileType','delimitedtext');
            %clusterTable=sortrows(clusterTable,'ch');
            qMetric=sortrows(qMetric,'maxChannels');
            %spike_templates = readNPY([pathToPhyResults filesep 'spike_templates.npy']);
            spikeTimes_samplesND = spike_times(~dspikes); %exclude duplicate spikes
            spikeTemplatesND = spike_clusters(~dspikes);%exclude duplicate spikes
            labelVec = {'noise','good','mua','non-somatic'};
            label = arrayfun(@(x) labelVec{x}, unitType+1, 'UniformOutput', false); %Replace values of label by their meaning
            %neuronAmp=clusterTable.amp;
               
            GoodUtemplate = find(unitType ~= 0); %Selects everything except noise units

            ic = zeros(5,numel(GoodUtemplate)); 

            t = cell(1,numel(GoodUtemplate));
            currentIdx=0;prevCh=-1;

            for i = 1:numel(GoodUtemplate)

                if ~isempty(spikeTimes_samplesND(spikeTemplatesND == GoodUtemplate(i)))

                    tm = spikeTimes_samplesND(spikeTemplatesND == GoodUtemplate(i));

                    t{i} = uint64(tm');

                    ic(1,i) = qMetric.maxChannels(GoodUtemplate(i),:);
                    ic(3,i)=currentIdx+1;
                    ic(4,i)=currentIdx+numel(t{i});
                    ic(5,i)=qMetric.phy_clusterID(GoodUtemplate(i),:);

                    if prevCh==ic(1,i)
                        ic(2,i)=ic(2,i-1)+1;
                    else
                        ic(2,i)=1;
                    end
                    prevCh=ic(1,i);
                    currentIdx=ic(4,i);

                end
            end

            nSpks=cellfun(@length,t(find(~cellfun(@isempty,t))));

            t=double(cell2mat(t(find(~cellfun(@isempty,t)))))/(obj.samplingFrequency(1)/1000);
            dic = ic(:,ic(1,:)~=0);

            fprintf('Saving results to %s\n',saveFileValid);
            save(saveFileAll,'t','ic','label','nSpks'); %save full spikes including noise

            %neuronAmp=neuronAmp(GoodUtemplate);
            
            save(saveFileValid,'t','ic','label','nSpks');
            
            if nargout==1 %if output is needed and calculation was needed (no saved file existing).
                spkData=load(saveFileValid);
            end
        end
        
        function []=convertLayoutKSort(obj,outputFile,badChannels)
            if nargin<2
                if iscell(obj.recordingDir)
                    outputFile=fullfile(obj.recordingDir{1}, 'chanMap.mat');
                else
                    outputFile=fullfile(obj.recordingDir, 'chanMap.mat');
                end
            end
            if nargin<3
                badChannels=[];
            end
            
            % here I know a priori what order my channels are in.  So I just manually
            % make a list of channel indices (and give an index to dead channels too). chanMap(1) is the row in the raw binary
            % file for the first channel. chanMap(1:2) = [33 34] in my case, which happen to be dead channels.
            
            chanMap = obj.channelNumbers;
            
            % the first thing Kilosort does is reorder the data with data = data(chanMap, :).
            % Now we declare which channels are "connected" in this normal ordering,
            % meaning not dead or used for non-ephys data
            
            badChannels=[];
            connected = true(numel(chanMap), 1);
            if ~isempty(badChannels)
                connected(badChannels) = false;
            end
            
            % now we define the horizontal (x) and vertical (y) coordinates of these
            % 34 channels. For dead or nonephys channels the values won't matter. Again
            % I will take this information from the specifications of the probe. These
            % are in um here, but the absolute scaling doesn't really matter in the
            % algorithm.
            
            xcoords = obj.chLayoutPositions(1,:);
            ycoords = obj.chLayoutPositions(2,:);
            
            % Often, multi-shank probes or tetrodes will be organized into groups of
            % channels that cannot possibly share spikes with the rest of the probe. This helps
            % the algorithm discard noisy templates shared across groups. In
            % this case, we set kcoords to indicate which group the channel belongs to.
            % In our case all channels are on the same shank in a single group so we
            % assign them all to group 1.
            
            kcoords = true(numel(chanMap), 1);
            
            % at this point in Kilosort we do data = data(connected, :), ycoords =
            % ycoords(connected), xcoords = xcoords(connected) and kcoords =
            % kcoords(connected) and no more channel map information is needed (in particular
            % no "adjacency graphs" like in KlustaKwik).
            % Now we can save our channel map for the eMouse.
            
            % would be good to also save the sampling frequency here
            fs = obj.samplingFrequency;
            
            save(outputFile, 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'fs');
            
            disp(['Channel information saved in : ' outputFile]);
            
        end

        function []=convertLayoutJRClust(obj,padSize,outputName)
            %convertLayouteJRClust(obj,padSize,outputName)
            %Make probe (.prb) file for using with jrclust
            %pad size - [height (y),widht (x)]
            %outputName - a name of the prb file (e.g. probe1.prb)        nCh=numel(~isnan(Enp(1,:)));
            %fid=fopen('layout_100_12x12.prb','w');
            if nargin<3
                outputName=[obj.recordingDir filesep obj.layoutName '_JRC.prb'];
            end
            if nargin<2
                error('Pad size must be entered as an external parameter');
            end
            if iscell(outputName) % temporary fix for outputname in multi file mode
                    Index = find(contains(outputName,'layout'));
                    temp = strsplit(outputName{Index}, '.mat');
                    new_outputname = strcat(temp{1},'_JRC.prb');
                    obj.layoutName = new_outputname;
                    fid=fopen(new_outputname,'w');
            else
                fid=fopen(outputName,'w');
            end
                        
            nCh=size(obj.chLayoutPositions,2);
            fprintf(fid, 'channels = [1:%d];\n\n',nCh);
            fprintf(fid, 'geometry = [%.1f,%.1f',obj.chLayoutPositions(1,1),obj.chLayoutPositions(2,1));
            for i=2:nCh
                fprintf(fid,';%.1f,%.1f',obj.chLayoutPositions(1,i),obj.chLayoutPositions(2,i));
            end
            fprintf(fid, '];\n\n');
            
            fprintf(fid, 'pad = [%.1f,%.1f];\n\n',padSize(1),padSize(2));
            
            fprintf(fid, 'cviShank = {1:%d};',nCh);
            fclose(fid);
        end
        
        function generateChannelMapFile(obj,electrodeName)
            overwrite=1;
            layoutFile=dir([obj.recordingDir filesep '*.chMap']);
            [layoutDir]=fileparts(which('layout_40_16x2_FlexLin.mat')); %an example of a file
            if ~isempty(layoutFile) && ~overwrite
                fprintf('\nLayout file already exists - %s\n',layoutFile.name);
                return;
            elseif ~exist('electrodeName','var')
                [layoutFile] = uigetfile([layoutDir filesep 'layout_*.mat'],'Select the electrode layout file');
                electrodeName=layoutFile(8:end-4);
            end
            fullFile=[layoutDir filesep 'layout_' electrodeName '.mat'];
            if isfile(fullFile)
                fid=fopen([obj.recordingDir filesep 'electrode.chMap'],'w');
                fprintf(fid,electrodeName);
                fclose(fid);
            else
                fprintf('Could not find %s!\n, check that electrode name was entered correctly or run without an electrode name and choose from the list.\n',fullFile);
            end
        end

        function convert2Binary(obj,targetFileBase,varargin)
            %Converts recording data to a simple binary format. 

            %parsing input variables
            parseObj = inputParser;
            addRequired(parseObj,'targetFileBase');
            addParameter(parseObj,'electrodeCh',obj.channelNumbers); %the electrode channel numbers to extract
            addParameter(parseObj,'analogCh',obj.analogChannelNumbers);%the analog channel numbers to extract
            addParameter(parseObj,'newElectrodeChNumbers',[]);%the new electrode channel numbers in the generated binary file
            addParameter(parseObj,'newAnalogChNumbers',[]);%the new analog channel numbers in the generated binary file
            addParameter(parseObj,'timeLimitsMs',[],@isnumeric);%[2 X N] - a vector with start and end times to extract [ms] - if N>2 combines together a few segments 
            addParameter(parseObj,'targetDataType','int16',@isstring); %class of electrode/analog data sample
            addParameter(parseObj,'overwriteBinaryData',0,@isnumeric); %if true overwrites data even if files exist
            addParameter(parseObj,'overwriteTriggerData',0,@isnumeric); %if true overwrites triggers even if files exist
            addParameter(parseObj,'chunkSize',2*60*1000,@isnumeric); %msec
            addParameter(parseObj,'downSample2LFP',0,@isnumeric); %if true downsamples by a factor of 250 to get LFP in 20KHz sampled signals
            addParameter(parseObj,'medianFilterGroup',[],@iscell);% {1 X N} - cell array of groups of electrode channel numbers for calculating joint median (if empty, does not filter)
            
            addParameter(parseObj,'inputParams',0,@isnumeric);
            parseObj.parse(targetFileBase,varargin{:});
            if parseObj.Results.inputParams %if true, plots all possible input params to the function
                disp(parseObj.Results);
                return;
            end
            par=parseObj.Results;  %make parameter structure

            %find orginal file recording name and directory
            %{
            if iscell(obj.recordingDir)
                recordingDir=obj.recordingDir{1};
                recordingName = strjoin(cellfun(@(x) x(1:end-3), obj.dataFileNames,'UniformOutput',false),'-');
            else
                recordingDir=obj.recordingDir;
                recordingName = obj.recordingName;
            end
            %}
            %Using default file name if an empty target is entered.
            if isempty(targetFileBase)
                if isempty(obj.recordingName)
                    targetFileBase=[obj.recordingDir filesep 'binaryRec.bin'];
                else
                    targetFileBase=[obj.recordingDir filesep obj.recordingName '.bin'];
                end
                fprintf('Target file name not entered! Using this: %s\n',targetFileBase);
            end
            %Verify that target extension is "bin" or "dat"
            if ~any(strcmp(targetFileBase(end-3:end),{'.dat','.bin'}))
                error('input file should have a ''.dat/.bin'' extension');
            end

            %create a folder if target folder does not exist
            [folderName,FileName]=fileparts(targetFileBase);
            if ~isfolder(folderName)
                fprintf('Notice the folder you entered does not exist!!! Creating new folder: %s\n',folderName);
                mkdir(folderName);
            end

            %determine conversion scheme based on the target file data class
            switch obj.datatype
                case par.targetDataType
                    zeroValue=0;
                    convertDataType=false; %switch from uint to int
                    convertTo16=false; %switch from 32int to 16int
                    useDouble=false; %other unknown type, us the data converted first to uV
                case 'uint16'
                    fprintf('Recording data type different from target data type!!!!\nConverting from %s to %s!\n',obj.datatype,par.targetDataType);
                    bits=16;
                    zeroValue=2^bits/2;
                    convertDataType=true;
                    convertTo16=false;
                    useDouble=false;
                case 'int32'
                    fprintf('Recording data type different from target data type!!!!\nConverting from %s to %s!\n',obj.datatype,par.targetDataType);
                    zeroValue=0;
                    convertTo16=true;
                    convertDataType=false;
                    useDouble=false;
                case 'uint32'
                    fprintf('Recording data type different from target data type!!!!\nConverting from %s to %s!\n',obj.datatype,par.targetDataType);
                    bits=32;
                    zeroValue=2^bits/2;
                    convertDataType=true;
                    convertTo16=true;
                    useDouble=false;
                otherwise
                    useDouble=true;
            end

            %Determine start and end of each processing chunk
            if isempty(par.timeLimitsMs)
                startTimes=0:par.chunkSize:obj.recordingDuration_ms;
                endTimes=[startTimes(2:end) obj.recordingDuration_ms];
            else
                %checks that the format of timeLimitsMs is correct and converts it if possible or gives an error.
                nTimeBlocks=size(par.timeLimitsMs);
                if nTimeBlocks(1)>2
                    error('The timeLimitsMs should be of size [2 x N] - where N is the number of chuncks!');
                elseif nTimeBlocks(1)==1
                    par.timeLimitsMs=par.timeLimitsMs';
                end
                %calculate start and end times for all sessions.
                startTimes=[];endTimes=[];
                for i=1:size(par.timeLimitsMs,2)
                    if par.timeLimitsMs(2,i)-(par.timeLimitsMs(1,i))>=par.chunkSize
                        tmpStart=par.timeLimitsMs(1,i):par.chunkSize:(par.timeLimitsMs(2,i)-par.chunkSize);
                    else
                        tmpStart=par.timeLimitsMs(1,i);
                    end
                    startTimes=[startTimes tmpStart];
                    endTimes=[endTimes tmpStart(2:end) par.timeLimitsMs(2,i)];
                end
                if ~issorted(startTimes) || ~issorted(endTimes)
                    error('timeLimitsMs must be sorted!!!');
                end
            end

            %check channel numbers
            if isempty(par.electrodeCh)
                electrodeCh=obj.channelNumbers;
            else
                electrodeCh=par.electrodeCh;
            end
            if isempty(par.analogCh)
                analogCh=obj.analogChannelNumbers;
            else
                analogCh=par.analogCh;
            end
            nAnalog=numel(analogCh);

            if ~isempty(par.newElectrodeChNumbers)
                if numel(par.newElectrodeChNames)~=numel(electrodeCh)
                    error('The number of new electrode channels and original electrode channels should be equal!!!');
                end
            else
                par.newElectrodeChNumbers=electrodeCh;
            end
            if ~isempty(par.newAnalogChNumbers)
                if numel(par.newAnalogChNumbers)~=numel(analogCh)
                    error('The number of new analog channels and original analogf= channels should be equal!!!');
                end
            else
                par.newAnalogChNumbers=analogCh;
            end

            %check if to plot all possible input params.
            for i=1:numel(par.medianFilterGroup)
                [~,pGroup{i}]=intersect(electrodeCh,par.medianFilterGroup{i});
            end

            if par.downSample2LFP
                %Creates decimation filter
                F=filterData(obj.samplingFrequency);
                F.downSamplingFactor=obj.samplingFrequency/250;
                F=F.designDownSample;
                F.padding=true;
                outSamplingFrequency=F.filteredSamplingFrequency;
            else
                outSamplingFrequency=obj.samplingFrequency(1);
            end
            %start extracting data
            tic;

            %Verifys that the target file does not exist to not overwrite files.
            if ~exist(targetFileBase,'file') || par.overwriteBinaryData
                try %convert electrode data to binary

                    %open data files
                    fid = fopen(targetFileBase, 'w+');
                    if nAnalog>0
                        fidA = fopen([targetFileBase(1:end-4) '_Analog.bin'], 'w+');
                    else
                        fidA=[];
                    end


                    %verify that dataRecording property is set such that all data in converted to double upon extraction
                    tempConvertData2Double=obj.convertData2Double; %loads data in double format
                    obj.convertData2Double=useDouble;

                    fprintf('\nConverting blocks to binary %s format(/%d) : ',par.targetDataType,numel(startTimes));
                    nDigits=0;
                    for j=1:numel(startTimes)
                        fprintf([repmat('\b',1,nDigits) '%d'],j);nDigits=length(num2str(j));
                        if ~useDouble
                            if ~convertDataType
                                data=squeeze(obj.getData(electrodeCh,startTimes(j),endTimes(j)-startTimes(j)));
                                if nAnalog>0
                                    dataAnalog=squeeze(obj.getAnalogData(analogCh,startTimes(j),endTimes(j)-startTimes(j)));
                                end
                                if convertTo16 %data is 32bit,convert to 16 with saturation
                                    saturation=32767; %2^15 -1 maximal value of 16bit
                                    data(data>saturation)=saturation;
                                    data(data<(-saturation))=-saturation;
                                    data=int16(data);
                                    if nAnalog>0
                                        dataAnalog(dataAnalog>saturation)=saturation;
                                        dataAnalog(dataAnalog<(-saturation))=-saturation;
                                        dataAnalog=int16(dataAnalog);
                                    end
                                end
                            else %convert from uint to int
                                if ~convertTo16
                                    data=int16(int32(squeeze(obj.getData(electrodeCh,startTimes(j),endTimes(j)-startTimes(j))))-zeroValue);
                                    if nAnalog>0
                                        dataAnalog=int16(int32(squeeze(obj.getAnalogData(analogCh,startTimes(j),endTimes(j)-startTimes(j))))-zeroValue);
                                    end
                                else
                                    %convert to signed
                                    saturation=32767; %2^15 -1 maximal value of 16bit
                                    data=int32(int64(squeeze(obj.getData(electrodeCh,startTimes(j),endTimes(j)-startTimes(j))))-zeroValue);
                                    %convert to 16 w saturation
                                    data(data>saturation)=saturation;
                                    data(data<(-saturation))=-saturation;
                                    data=int16(data);

                                    if nAnalog>0
                                        dataAnalog=int32(int64(squeeze(obj.getAnalogData(analogCh,startTimes(j),endTimes(j)-startTimes(j))))-zeroValue);
                                        dataAnalog(dataAnalog>saturation)=saturation;
                                        dataAnalog(dataAnalog<(-saturation))=-saturation;
                                        dataAnalog=int16(dataAnalog);
                                    end
                                end
                            end
                        else %unknown data type - get data in uV and convert to 16bit
                            saturation=32767; %2^15 -1 maximal value of 16bit
                            data=squeeze(obj.getData(electrodeCh,startTimes(j),endTimes(j)-startTimes(j)));
                            data=(data+obj.ZeroADValue)/obj.MicrovoltsPerAD(1);
                            data(data>saturation)=saturation;
                            data(data<(-saturation))=-saturation;
                            data=int16(data);
                            if nAnalog>0
                                dataAnalog=squeeze(obj.getAnalogData(analogCh,startTimes(j),endTimes(j)-startTimes(j)));
                                dataAnalog=(dataAnalog+obj.ZeroADValue)/obj.MicrovoltsPerADAnalog(1);
                                dataAnalog(dataAnalog>saturation)=saturation;
                                dataAnalog(dataAnalog<(-saturation))=-saturation;
                                dataAnalog=int16(dataAnalog);
                            end
                        end
                        %median filters are performed only on electrode channels (not on analog channels)
                        if ~isempty(par.medianFilterGroup)
                            for i=1:numel(par.medianFilterGroup)
                                data(pGroup{i},:)=bsxfun(@minus,data(pGroup{i},:),median(data(pGroup{i},:)));
                            end
                        end
                        pause(0.0001);
                        if par.downSample2LFP
                            data=int16(shiftdim(F.getFilteredData(double(shiftdim(data,-1))),1));
                        end
                        fwrite(fid, data, ['*' par.targetDataType]);
                        if nAnalog>0
                            fwrite(fidA, dataAnalog, ['*' par.targetDataType]);
                        end
                    end
                    fclose(fid);
                    if ~isempty(fidA)
                        fclose(fidA);
                    end
                    fprintf('\nConversion complete!\n');
                    obj.convertData2Double=tempConvertData2Double; %return value to what it was

                catch ME
                    fclose(fid);
                    if ~isempty(fidA)
                        fclose(fidA);
                    end
                    obj.convertData2Double=tempConvertData2Double; %return value to what it was
                    fprintf('Error in conversion! - closing files...\n');
                    rethrow(ME)
                end
            else
                disp('The chosed binary data file name already exists!!!! Skipping this file!');
            end

            %extracts trigger information
            fidD=[];
            triggerFile=[targetFileBase(1:end-4) '_Triggers.bin'];
            if ~exist(triggerFile,'file') || par.overwriteTriggerData
                try %convert electrode data to binary
                    fprintf('\nConverting digital trigger file...\n');
                    T=obj.getTrigger;
                    nT=cellfun(@(x) numel(x),T);
                    pT=find(nT>0);
                    for k=1:numel(pT)
                        if size(T{pT(k)},1)>1
                            T{pT(k)}=T{pT(k)}';
                            fprintf('Notice! The input triggers are not all [1 x N] vectors. Changing and continuing...\n')
                        end
                    end
                    %correct time stamps according to cuts
                    TNew=cell(1,numel(nT));
                    accumDelay=cumsum([startTimes(1) startTimes(2:end)-endTimes(1:end-1)]);
                    for j=1:numel(startTimes)
                        for k=1:numel(pT)
                            pTmp=find(T{pT(k)}>=startTimes(j) & T{pT(k)}<endTimes(j));
                            TNew{pT(k)}=[TNew{pT(k)} (T{pT(k)}(pTmp(:))-accumDelay(j))*obj.samplingFrequency(1)/1000];
                        end
                    end
                    nT=cellfun(@(x) numel(x),TNew);

                    fidD = fopen(triggerFile, 'w');
                    fwrite(fidD,uint32(nT+1),'*uint32');%write the number of triggers from each type.
                    for j=1:numel(pT)
                        fwrite(fidD, uint32(TNew{pT(j)}+1),'*uint32');
                    end
                    fclose(fidD);
                catch ME
                    if ~isempty(fidD)
                        fclose(fidD);
                    end
                    nT = [];
                    disp('No triggers found! Trigger file not created.\n');
                    rethrow(ME)
                end
            else
                disp('The chosen trigger file name already exists!!!! Skipping this file!');
            end

            %generates meta data file
            metaDataFile=[targetFileBase(1:end-4) '_meta.txt'];
            if ~exist(metaDataFile,'file')
                fidM=fopen(metaDataFile,'w');
                fprintf(fidM,'nSavedChans = %d\n',numel(electrodeCh));
                fprintf(fidM,'sRateHz = %.12f\n',outSamplingFrequency);
                fprintf(fidM,'nChans = %d\n',numel(electrodeCh));
                outputstr = ['%d' repmat(',%d', 1, numel(par.newElectrodeChNumbers)-1)]; % replicate it to match the number of columns
                fprintf(fidM,['channelNumbers = ', outputstr, '\n'], par.newElectrodeChNumbers);

                fprintf(fidM,'nAnalogChans = %d\n',numel(obj.analogChannelNumbers));
                if nAnalog>0
                    fprintf(fidM,'sRateAnalogHz = %.12f\n',obj.samplingFrequencyAnalog(1));
                    outputstr = ['%d' repmat(',%d', 1, numel(par.newAnalogChNumbers)-1)]; % replicate it to match the number of columns
                    fprintf(fidM,['channelNumbersAnalog = ', outputstr, '\n'], par.newAnalogChNumbers);
                    fprintf(fidM,'scaleAnalog = %.12f\n',obj.MicrovoltsPerADAnalog(1));
                    fprintf(fidM,'zeroADValueAnalog = %.12f\n',obj.ZeroADValueAnalog(1));
                else
                    fprintf(fidM,'sRateAnalogHz = []\n');
                    fprintf(fidM,'channelNumbersAnalog = []\n');
                    fprintf(fidM,'scaleAnalog = []\n');
                    fprintf(fidM,'zeroADValueAnalog = []\n');
                end

                fprintf(fidM,'nTriggerChans = %d\n',numel(nT));
                fprintf(fidM,'vcDataType = %s\n',par.targetDataType);

                fprintf(fidM,'scale = %.12f\n',obj.MicrovoltsPerAD(1));
                fprintf(fidM,'zeroADValue = %.12f\n',obj.ZeroADValue(1));


                fprintf(fidM,'vcProbe = %s\n',obj.layoutName);
                fclose(fidM);
            else
                disp(['Meta data file: ' metaDataFile ' alreday exists!, please first delete']);
            end
            toc;
        end
        
        function obj=getRecordingFiles(obj,recordingFile,fileExtension)
            %Get directory with data files
            %Usage: obj = getRecordingFiles(obj,recordingFile,fileExtension)
            %if no recording file is entered lauches a GUI
            %if no file extension entered, a directory is chosen rather than a specific files (for example for neuralynx recordings)
            
            %If no files were entered open GUI for choosing a file or a directory else get the files entered
            if ~isempty(recordingFile) %if directory with data was not entered open get directory GUI
                
                obj.multifileMode=iscell(recordingFile);
                if obj.multifileMode,singleRecordingFile=recordingFile{1};else singleRecordingFile=recordingFile;end
                if isdir(singleRecordingFile)
                    obj.folderMode=true; %a folder is chosen and the files inside examined
                else
                    obj.folderMode=false; %a file or list of files is selected
                end
                
                if ~obj.folderMode
                    if ~obj.multifileMode
                        recordingFile={recordingFile};
                    end
                    obj.nRecordings=numel(recordingFile);
                    for i=1:obj.nRecordings
                        [pathstr{i}, name{i}, ext] = fileparts(recordingFile{i});
                        obj.dataFileNames{i}=[name{i} ext];
                        if ~exist([pathstr{i} filesep obj.dataFileNames{i}],'file')
                            disp(['Searching for recording file: ' [pathstr{i} filesep obj.dataFileNames{i}]]);
                            error('Object was not constructed since no valid recording file name was chosen');
                        end
                    end
                else
                    if ~obj.multifileMode
                        [pathstr, name] = fileparts(recordingFile);
                        obj.dataFileNames{1}=recordingFile;
                    else
                        [pathstr, name] = cellfun(@(x) fileparts(x),recordingFile,'UniformOutput',0);
                        obj.dataFileNames=recordingFile;
                    end
                end
                
                if isempty(pathstr) %in case the file is in the current directory
                    if ispc
                        obj.recordingDir=[cd filesep];
                    end
                else
                    if iscell(pathstr) && all(strcmp(pathstr,pathstr{1}))
                        pathstr=pathstr{1};
                    end
                    if iscell(name) && all(strcmp(name,name{1}))
                        name=name{1};
                    end
                    obj.recordingDir=pathstr;
                    if ispc
                        if ~iscell(obj.recordingDir)
                            obj.recordingDir=[obj.recordingDir filesep];
                        else
                            obj.recordingDir=cellfun(@(x) [x filesep],obj.recordingDir,'UniformOutput',0);
                        end
                    end
                end
                
                if iscell(obj.recordingDir) %some of the condition below can be removed
                    if ~isdir(obj.recordingDir{1})
                        error('Object was not constructed since no valid folder was choosen');
                    end
                else
                    if ~isdir(obj.recordingDir)
                        error('Object was not constructed since no valid folder of file name were choosen!!!');
                    end
                end
            else %if directory with data was not entered open get directory GUI
                if ~obj.folderMode
                    [obj.dataFileNames,obj.recordingDir]= uigetfile(['*.' fileExtension],['Choose the ' fileExtension ' file'],obj.defaultLocalDir,'MultiSelect','on');
                    if ~iscell(obj.dataFileNames)
                        obj.dataFileNames={obj.dataFileNames};
                    end
                    if obj.dataFileNames{1}==0 %no folder chosen
                        disp('Object was not constructed since no folder was choosen');
                        return;
                    end
                    obj.nRecordings=numel(obj.dataFileNames);
                    if obj.nRecordings>1
                        obj.multifileMode=true; 
                    end
                else
                    [obj.recordingDir]= uigetdir(obj.defaultLocalDir,'Choose the data folder');
                    [pathstr, name] = fileparts(obj.recordingDir);
                    obj.recordingDir=[pathstr filesep];
                    obj.multifileMode=false;
                end
            end
            if ~obj.folderMode
                [~,obj.recordingName]=fileparts(obj.dataFileNames{1});
            else
                obj.recordingName=name;
            end
            %for some class like OERecording a speicific metadata file name is give, if not use this name for metadata
            if isempty(obj.metaDataFile)
                if ~iscell(obj.recordingDir)
                    obj.metaDataFile=[obj.recordingDir filesep obj.recordingName '_metaData.mat'];
                else
                    obj.metaDataFile=[obj.recordingDir{1} filesep obj.recordingName '_metaData.mat'];
                end
            end
            if obj.recordingDir(end)=='/'
                obj.recordingDir=obj.recordingDir(1:end-1);
            end
        end

    end
end