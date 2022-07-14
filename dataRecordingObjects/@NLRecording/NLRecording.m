classdef NLRecording < dataRecording
    properties
        
        recordSize
        ADMaxValue
        ADBitVolts
        inputRange
        inputInverted
         
        triggerFileName %The name (including full path) of the trigger file in Neuralynx recording
        recordingStartTime_us
        
        recordLength_us
        sample_us
        nRecords
        maxTTLBit = 9;
        
        allTimeStamps_us
        recFileOpen
        recFileClosed
        
        isunix
    end
    
    properties (SetAccess=protected) %these properties are not considered as meta data and have to be loaded each time
        fid %file identifiers for loading data
        channelFiles
    end
    
    properties (Constant)
        defaultLocalDir='C:\Users\Tulip\Documents\Academic\Post-Doc\Experiments'; %Default directory from which search starts
        headerSizeByte=16384;
        fileExtension='ncs';
        dataSamplesPerRecord=512;
        bytesPerRec=1044; %64+32+32+32+512*16 bits = 1044 bytes
        signalBits = 16; %the quantization of the sampling card
        numberOfCharFromEndToBaseName=7;
    end
    
    methods
        
        function [V_uV,t_ms]=getData(obj,channels,startTime_ms,window_ms)
            %Extract Neuralynx recording data from file to memory
            %Usage: [V_uV,t_ms]=obj.getData(channels,startTime_ms,window_ms);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        t_ms - A time vector relative to recording start (t=0 at start)
            startTime_us=startTime_ms*1e3;
            windowSamples=round(window_ms/1e3*obj.samplingFrequency(1));
            window_us=window_ms*1e3;
            nWindows=numel(startTime_us);
            
            if isempty(channels) %if no channels are entered, get all channels
                channels=obj.channelNumbers;
            end
            nCh=numel(channels);
            
            V_uV=zeros(windowSamples,nWindows,nCh,'int16'); %initialize waveform matrix
            
            %generate time stamps for block waveform extraction and extract waveforms from file
            
            for i=1:nWindows
                pSingleTrialTimeStamps{i}=find(obj.allTimeStamps_us>=startTime_us(i)-obj.recordLength_us & obj.allTimeStamps_us<(startTime_us(i)+window_us));
                singleTrialTimeStamps=obj.allTimeStamps_us(pSingleTrialTimeStamps{i});
                recordsPerTrial(i)=numel(singleTrialTimeStamps);
                timeIdx=bsxfun(@plus,(1:obj.dataSamplesPerRecord)*obj.sample_us(1),singleTrialTimeStamps);
                pRecIdx{i,:}=timeIdx>=startTime_us(i) & timeIdx<(startTime_us(i)+window_us);
                timeIdx=timeIdx';
                pOutIdx{i,1}=floor((timeIdx(pRecIdx{i,:}')-startTime_us(i))./obj.sample_us(1)+1)+windowSamples*(i-1);
            end
            pRecIdx=cell2mat(pRecIdx);
            pOutIdx=cell2mat(pOutIdx);
            
            for i=1:nCh
                data=zeros(size(pRecIdx),'int16')';currRec=1;
                for j=1:nWindows
                    if ~isempty(pSingleTrialTimeStamps{j})
                        %obj.fid(obj.n2s(channels(i)))=fopen([obj.recordingDir filesep obj.channelFiles{pCh(i)}],'r');
                        fseek(obj.fid(obj.n2s(channels(i))),obj.headerSizeByte+(pSingleTrialTimeStamps{j}(1)-1)*obj.bytesPerRec+20,'bof'); %(64+32+32+32)/8=20
                        data(:,(1:recordsPerTrial(j))+currRec-1)=fread(obj.fid(obj.n2s(channels(i))), [obj.dataSamplesPerRecord recordsPerTrial(j)], '512*int16',obj.bytesPerRec-obj.dataSamplesPerRecord*16/8);
                        currRec=currRec+recordsPerTrial(j);
                    else
                        disp('requested time stamp outside recording range!');
                    end
                end
                V_uV(pOutIdx+(i-1)*nWindows*windowSamples)=data(pRecIdx');
            end
            
            if obj.convertData2Double
                V_uV = permute((2*(0.5-obj.inputInverted(1)))*(double(V_uV) - obj.ZeroADValue) * obj.MicrovoltsPerAD ,[3 2 1]);
            else
                V_uV = permute(V_uV,[3 2 1]);
            end
            
            if nargout==2
                t_ms=(1:windowSamples)*(1e3/obj.samplingFrequency(1));
            end
        end
        
        function [T_ms,chNumber,chName]=getTrigger(obj,startTime_ms,window_ms,name)
            %Extract triggers from file Neuralynx recording
            %Usage : [T_ms]=obj.getTrigger(name,startTime_ms,endTime_ms)
            %Input : name - which bit to extract for time stamps (out of 8,default = first bit, 1)
            %        startTime_ms - start time [ms].
            %        window_ms - the window duration [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %
            %Output: T_ms - trigger times [ms] - different triggers are arranged in a cell array
            
            %read Data from *.nev files (NeuroLynx) - files containing the data
            %fieldSelectionTrigger: (1): Timestamps, (2): Event IDs, (3): TTLs, (4): Extras, (5): Event Strings
            %Bits are coded as decimal integers
            if nargin>=3
                if ~obj.isunix
                [timeStampsTrig, TTLsTrig,eventString] = Nlx2MatEV...
                    ([obj.recordingDir filesep obj.triggerFileName],[1 0 1 0 1],0,4,1e3*[startTime_ms; startTime_ms+window_ms]+obj.recordingStartTime_us);
                else
                [timeStampsTrig, TTLsTrig,eventString] = Nlx2MatEV_v3...
                    ([obj.recordingDir filesep obj.triggerFileName],[1 0 1 0 1],0,4,1e3*[startTime_ms; startTime_ms+window_ms]+obj.recordingStartTime_us);
                end
            else
                if ~obj.isunix
                [timeStampsTrig, TTLsTrig,eventString] = Nlx2MatEV...
                    ([obj.recordingDir filesep obj.triggerFileName],[1 0 1 0 1],0,1,[]);
                else
                [timeStampsTrig, TTLsTrig,eventString] = Nlx2MatEV_v3...
                    ([obj.recordingDir filesep obj.triggerFileName],[1 0 1 0 1],0,1,[]);
                end
            end
            
            %calculate bit from integer sequence
            maxTTLBit=8; %the max bit for checking regular TTLs, bits above this value will be checked for message TTLs
            
            pSimpleTTLBinary=TTLsTrig<=2^(maxTTLBit);
            pSimpleTTL=find(pSimpleTTLBinary);
            TTLBinaryCh=logical(fliplr(rem(floor(TTLsTrig(pSimpleTTL)'*pow2(1-maxTTLBit:0)),2)));
            
            %check if specific bits were selected for extraction
            if nargin~=4
                bits2Check=1:maxTTLBit;
            else
                bits2Check=name;
            end
            
            %extract bits time stamps and names
            for i=bits2Check%(sum(TTLBinaryCh)>=1)
                pTmp=pSimpleTTL(1+find(diff(TTLBinaryCh(:,i))~=0));
                T_ms{i}=(   timeStampsTrig(pTmp) - obj.recordingStartTime_us   )/1e3;
                chName{i}=eventString(pTmp);
                chNumber(i)=i;
            end
            
            %extract TTLs that were coded as integers and do not necessarily have an off signal
            if isempty(bits2Check)
                bits2Check=0;
            end
            complexTTLStart=max(bits2Check)+1;
            complexTTLs=unique(TTLsTrig(~pSimpleTTLBinary));
            for i=1:numel(complexTTLs)
                pTmp=(find(TTLsTrig==complexTTLs(i)));
                T_ms{i+complexTTLStart}=(timeStampsTrig(pTmp)-obj.recordingStartTime_us)/1e3;
                chName{i+complexTTLStart}=eventString(pTmp);
                chNumber(i+complexTTLStart)=complexTTLs(i);
            end
            
        end
        
        function obj=getStartRecordingTime(obj)
            fileID = fopen([obj.recordingDir filesep 'CheetahLogFile.txt']);
            logFile=textscan(fileID,'%s','Delimiter','-');
            fclose(fileID);
            pStart=find(strcmp(logFile{1},'Recording Started'));
            recordingStart_us=logFile{1}(pStart-2);
            timeStart=logFile{1}(pStart-3);
            obj.startDate=datestr(datenum(datestr(timeStart{1},'HH:MM:SS:FFF'),'HH:MM:SS')  +  datenum(0,0,0,0,0, (obj.recordingStartTime_us-str2num(recordingStart_us{1}))/1e6 )  ,'HH:MM:SS');
        end
        
        function obj=getHeaderInformation(obj) %get parameters from file header
            disp('Extractign header information');
            obj=getStartRecordingTime(obj); %get sync to computer clock from log file in start date file
            
            %get header and first time stamp from recordings
            obj.startDate=zeros(1,numel(obj.channelNumbers));
            obj.endDate=zeros(1,numel(obj.channelNumbers));
            for i=1:numel(obj.channelNumbers)
                fseek(obj.fid(i),0,'bof');
                hd=fread(obj.fid(i),obj.headerSizeByte);
                jumpline=[0;find(hd==13)];
                headerData{i}=cell(numel(jumpline)-1,1);
                for l=1:numel(jumpline)-1
                    str=char(hd(jumpline(l)+1:jumpline(l+1)-1)');
                    if numel(str)>1
                        headerData{i}{l}=char(hd(jumpline(l)+2:jumpline(l+1)-1)');
                    else
                        headerData{i}{l}='';
                    end
                    
                    strParts=strsplit(headerData{i}{l},' ');
                    switch(strParts{1})
                        case '##'
                            if strParts{3}(1)=='O' %'Opened'
                                obj.recFileOpen(i)=datenum([strParts{5} '-' strParts{7}],'mm/dd/yyyy-HH:MM:SS.FFF');
                            elseif strParts{3}(1)=='C' %'Closed'
                                if strParts{4}(1)~='F'
                                    obj.recFileClosed(i)=datenum([strParts{5} '-' strParts{7}],'mm/dd/yyyy-HH:MM:SS.FFF');
                                else
                                    obj.recFileClosed(i)=Inf;
                                end
                            end
                        case '-RecordSize'
                            obj.recordSize(i)=str2num(strParts{2});
                        case '-AcqEntName'
                            obj.channelNames{i}=strParts{2};
                        case '-SamplingFrequency'
                            obj.samplingFrequency(i)=str2num(strParts{2});
                            obj.sample_us(i)=1/obj.samplingFrequency(i)*1e6;
                        case '-ADMaxValue'
                            obj.ADMaxValue(i)=str2num(strParts{2});
                        case '-ADBitVolts'
                            obj.ADBitVolts(i)=str2num(strParts{2});
                        case '-InputRange'
                            obj.inputRange(i)=str2num(strParts{2});
                        case '-InputInverted'
                            if strParts{2}(1)=='F'
                                obj.inputInverted(i)=0;
                            else
                                obj.inputInverted(i)=1;
                            end
                        case '-DspLowCutFrequency'
                            obj.dspLowCutFrequency(i)=str2num(strParts{2});
                        case '-DspHighCutFrequency'
                            obj.dspHighCutFrequency(i)=str2num(strParts{2});
                    end
                end
            end
            obj.ZeroADValue=0;
            obj.MicrovoltsPerAD=obj.ADBitVolts(1)*1e6;
        end
        
        function obj=closeOpenFiles(obj) %clear all open file handles
            for i=1:numel(obj.fid)
                fclose(obj.fid(i));
            end
        end
        
        function obj=extractMetaData(obj)
            obj=getHeaderInformation(obj);
            
            % the file is a 16K header with a bunch of CSC records afterwards.
            % each CSC record is (timestamp-ch-samplingRate-nSamples-data) 64+32+32+32+512*16 bits = 1044 bytes
            
            % go to end of file to get filesize check of the header info and the file size (I assume that all channel files have the same number of records and time stamps)
            disp('Extracting time stamp information from Neuralynx recording...');
            
            fseek(obj.fid(1), 0, 1);
            obj.nRecords=(ftell(obj.fid(1))-obj.headerSizeByte)/obj.recordSize(1);
            obj.recordLength_us=obj.dataSamplesPerRecord*obj.sample_us(1);
            
            fseek(obj.fid(1), obj.headerSizeByte, 'bof');
            obj.allTimeStamps_us=fread(obj.fid(1), obj.nRecords, 'int64',obj.recordSize(1)-64/8);
            obj.recordingStartTime_us=obj.allTimeStamps_us(1);
            obj.allTimeStamps_us=obj.allTimeStamps_us-obj.allTimeStamps_us(1);
            obj.recordingDuration_ms=obj.allTimeStamps_us(end)/1000;
            
            %load layout
            obj=obj.loadChLayout;
            if isempty(obj.chLayoutNumbers)
                disp('No layout file in directory, to add layout please add a *.chMap file to the data directory with the layout name');
            end
            
            triggerNames=dir([obj.recordingDir filesep '*.nev']);
            if numel(triggerNames)~=1
                disp('triggers not extracted - more than one trigger file in directory');
            else
                obj.triggerFileName=triggerNames.name;
            end
            disp('saving meta data');
            obj.saveMetaData;
        end
        
         function obj=getRecordingFiles(obj,recordingFile,fileExtension)
            %Get directory with data files
            %Usage: obj = getRecordingFiles(obj,recordingFile,fileExtension)
            %if no recording file is entered lauches a GUI
            %if no file extension entered, a directory is chosen rather than a specific files (for example for neuralynx recordings)
            
            %If no files were entered open GUI for choosing a file or a directory else get the files entered
            obj.multifileMode=false; %the ability to process multiple folders together was not implemented
            
            if isempty(recordingFile) %if directory with data was not entered open get directory GUI
                [recordingFile]= uigetdir(obj.defaultLocalDir,'Choose the data folder');
            end
            if iscell(recordingFile)
                if numel(recordingFile)>1
                    error('NLRecoding currently only works with one directory (multiFileMode not implemented)')
                else
                    recordingFile=recordingFile{1};
                end
            end
            if recordingFile(end)==filesep %required for separating recordingFile to pathstr and name (and not only pathstr with empty name)
                recordingFile(end)=[];
            end
            [pathstr, name] = fileparts(recordingFile);
            obj.recordingDir=recordingFile;
            obj.recordingName=name;
            obj.metaDataFile=[obj.recordingDir filesep obj.recordingName '_metaData'];
         end
         
    end
    
   

    
    methods (Hidden)
        %class constructor
        function obj = NLRecording(recordingFile)
            %Usage: obj = NeuraLynxRecording(recordingFile)

            %get data files
            if nargin==0
                recordingFile=[];
            elseif nargin>1
                disp('Object was not constructed since too many parameters were given at construction');
                return;
            end
            obj.isunix=isunix;
            
            obj.datatype='int16';
            obj.folderMode=true; %neuralynx always work in folder mode - a recording is defined by a folder
            obj=obj.getRecordingFiles(recordingFile);
            obj.channelFiles=dir([obj.recordingDir filesep '*.' obj.fileExtension]);
            obj.channelFiles={obj.channelFiles.name};
            
            %get channel information
            obj.channelNumbers=cell2mat(cellfun(@(x) str2double(x(4:end-4)),obj.channelFiles,'UniformOutput',0));
            pCh=isnan(obj.channelNumbers); %in some cases there are dummy files in the folder with ncs ending - these should be removed
            obj.channelNumbers(pCh)=[];
            obj.channelFiles(pCh)=[];
            
            [obj.channelNumbers,pCh]=sort(obj.channelNumbers);
            obj.channelFiles=obj.channelFiles(pCh);
            obj.n2s(obj.channelNumbers)=1:numel(obj.channelNumbers);
            
            %get file indentifiers for all channel files
            for i=1:numel(obj.channelFiles)
                obj.fid(i)=fopen([obj.recordingDir filesep obj.channelFiles{i}],'r');
            end
            
            if exist([obj.metaDataFile '.mat'],'file') && ~obj.overwriteMetaData
                obj=loadMetaData(obj);
            else
                obj=extractMetaData(obj);
            end
            
            obj=obj.loadChLayout;
            
        end
        
    end
    
end