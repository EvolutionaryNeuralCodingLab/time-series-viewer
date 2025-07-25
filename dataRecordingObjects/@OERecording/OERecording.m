classdef OERecording < dataRecording
    properties
        channelFiles
        channelFilesAnalog
        eventFiles
        evntFileSize
        fileSize
        fileSizeAnalog
        channelNameFromHeader
        channelTypeFromHeader
        channelNameFromHeaderAnalog
        channelTypeFromHeaderAnalog

        n2sA
        fileHeaders
        fileHeadersAnalog
        ADBitVolts
        recordSize
        
        sample_ms
        bufferSize;
        bufferSizeAnalog
        blockLength;
        blockLengthAnalog;
        
        startDateA
        blkCont
        nRecordsCont
        blkBytesCont
        dataDescriptionCont
        dataDescriptionContAnalog
        bytesPerRecCont
        recordLength
        
        blkEvnt
        nRecordsEvnt
        dataDescriptionEvnt
        bytesPerRecEvnt
        blkBytesEvnt
        eventFileName = 'all_channels.events'
        
        softwareVersion
        softwareVersionAnalog
        openEphyXMLData
        openEphyXMLStructureData
        
        allTimeStamps
        globalStartTime_ms
    end
    
    properties (SetAccess=protected) %these properties are not considered as meta data and have to be loaded each time
        fid %file identifiers for loading data
        fidA
        fidEvnt
    end
    
    properties (Constant)
        defaultLocalDir='C:\Users\Tulip\Documents\Academic\Post-Doc\Experiments'; %Default directory from which search starts
        headerSizeByte=1024;
        fileExtension='continuous';
        eventFileExtension='events';
        signalBits = 16; %the quantization of the sampling card
        dataSamplesPerRecord=1024;
        maxTTLBit = 9;
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
            windowSamples=round(window_ms/obj.sample_ms);
            nWindows=numel(startTime_ms);
            startTime_ms=round(startTime_ms/obj.sample_ms)*obj.sample_ms;
            window_ms=windowSamples*obj.sample_ms;
            
            if isempty(channels) %if no channels are entered, get all channels
                channels=obj.channelNumbers;
            end
            if ~all(ismember(channels,obj.channelNumbers))
                error('The entered channel number does not exist in the recording!');
            end
            nCh=numel(channels);
            
            V_uV=zeros(windowSamples,nWindows,nCh,obj.blkCont(4).Types); %initialize waveform matrix
            pOutIdx=cell(nWindows,1);
            pRecIdx=cell(nWindows,1);
            %generate time stamps for block waveform extraction and extract waveforms from file
            %clear pOutIdx
            for i=1:nWindows
                pSingleTrialTimeStamps{i}=find(obj.allTimeStamps>=startTime_ms(i)-obj.recordLength & obj.allTimeStamps<(startTime_ms(i)+window_ms)); %find relevant blocks in block list
                singleTrialTimeStamps=round(obj.allTimeStamps(pSingleTrialTimeStamps{i})/obj.sample_ms)*obj.sample_ms; %calculate time stamps in milliseconds
                recordsPerTrial(i)=numel(singleTrialTimeStamps);
                timeIdx=bsxfun(@plus,(0:obj.dataSamplesPerRecord-1)*obj.sample_ms,singleTrialTimeStamps); % create a matrix for the times of every sample in ms
                pRecIdx{i,:}=(timeIdx>=startTime_ms(i)) & timeIdx<(startTime_ms(i)+window_ms); %find time indices within the requested time window (chuncks are 1024 in size so they are usually cut for most time windows)
                
                %maybe better to replace with "if pOutIdx{i,1}==0"
                if sum(sum(pRecIdx{i,:}(:)))==windowSamples+1 %due to rounding issues, there may be an error when there is one sample too much - in this case the last sample is removed
                    pRecIdx{i,:}(1,find(pRecIdx{i,:}(1,:)==1,1,'first'))=false;
                end
                
                timeIdx=timeIdx';
                pOutIdx{i,1}=round((timeIdx(pRecIdx{i,:}')-startTime_ms(i))/obj.sample_ms)+windowSamples*(i-1); %round should not be changed to floor or ceil - it creates a weird artifact
            end
            %this solved a special case that may not be needed anymore - check if the future if can be removed
            if pOutIdx{1}(1)==0
                pOutIdx{1}=pOutIdx{1}+1;
            end
            pRecIdx=cell2mat(pRecIdx);
            pOutIdx=cell2mat(pOutIdx);
            
            %}
            
            for i=1:nCh
                data=zeros(size(pRecIdx),'int16')';currRec=1;
                for j=1:nWindows
                    if ~isempty(pSingleTrialTimeStamps{j})
                        fseek(obj.fid(obj.n2s(channels(i))),obj.headerSizeByte+(pSingleTrialTimeStamps{j}(1)-1)*obj.bytesPerRecCont+sum(obj.blkBytesCont(1:3)),'bof'); %(64+32+32+32)/8=20
                        data(:,(1:recordsPerTrial(j))+currRec-1)=fread(obj.fid(obj.n2s(channels(i))), [obj.dataSamplesPerRecord recordsPerTrial(j)], '1024*int16',obj.bytesPerRecCont - obj.blkBytesCont(4),'b');
                        currRec=currRec+recordsPerTrial(j);
                    else
                        disp('requested time stamp outside recording range!');
                    end
                end
                V_uV(pOutIdx+(i-1)*nWindows*windowSamples)=data(pRecIdx');
            end
            
            if obj.convertData2Double
                V_uV = permute(double(V_uV) * obj.MicrovoltsPerAD(1) ,[3 2 1]);
            else
                V_uV = permute(V_uV,[3 2 1]);
            end
            
            if nargout==2
                t_ms=(1:windowSamples)*(1e3/obj.samplingFrequency(1));
            end
        end
        
        function [V_uV,t_ms]=getAnalogData(obj,channels,startTime_ms,window_ms)
            %Extract Neuralynx recording data from file to memory
            %Usage: [V_uV,t_ms]=obj.getData(channels,startTime_ms,window_ms);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        t_ms - A time vector relative to recording start (t=0 at start)
            windowSamples=round(window_ms/obj.sample_ms);
            nWindows=numel(startTime_ms);
            startTime_ms=round(startTime_ms/obj.sample_ms)*obj.sample_ms;
            window_ms=windowSamples*obj.sample_ms;
            
            if isempty(channels) %if no channels are entered, get all channels
                channels=obj.analogChannelNumbers;
            end
            nCh=numel(channels);
            
            V_uV=zeros(windowSamples,nWindows,nCh,obj.blkCont(4).Types); %initialize waveform matrix
            pOutIdx=cell(nWindows,1);
            pRecIdx=cell(nWindows,1);
            %generate time stamps for block waveform extraction and extract waveforms from file
            %clear pOutIdx
            for i=1:nWindows
                pSingleTrialTimeStamps{i}=find(obj.allTimeStamps>=startTime_ms(i)-obj.recordLength & obj.allTimeStamps<(startTime_ms(i)+window_ms)); %find relevant blocks in block list
                singleTrialTimeStamps=round(obj.allTimeStamps(pSingleTrialTimeStamps{i})/obj.sample_ms)*obj.sample_ms; %calculate time stamps in milliseconds
                recordsPerTrial(i)=numel(singleTrialTimeStamps);
                timeIdx=bsxfun(@plus,(0:obj.dataSamplesPerRecord-1)*obj.sample_ms,singleTrialTimeStamps); % create a matrix for the times of every sample in ms
                pRecIdx{i,:}=(timeIdx>=startTime_ms(i)) & timeIdx<(startTime_ms(i)+window_ms); %find time indices within the requested time window (chuncks are 1024 in size so they are usually cut for most time windows)
                
                %maybe better to replace with "if pOutIdx{i,1}==0"
                if sum(sum(pRecIdx{i,:}(:)))==windowSamples+1 %due to rounding issues, there may be an error when there is one sample too much - in this case the last sample is removed
                    pRecIdx{i,:}(1,find(pRecIdx{i,:}(1,:)==1,1,'first'))=false;
                end
                
                timeIdx=timeIdx';
                pOutIdx{i,1}=round((timeIdx(pRecIdx{i,:}')-startTime_ms(i))/obj.sample_ms)+windowSamples*(i-1); %round should not be changed to floor or ceil - it creates a weird artifact
            end
            %this solved a special case that may not be needed anymore - check if the future if can be removed
            if pOutIdx{1}(1)==0
                pOutIdx{1}=pOutIdx{1}+1;
            end
            pRecIdx=cell2mat(pRecIdx);
            pOutIdx=cell2mat(pOutIdx);
            
            %{
            if pOutIdx(1)==0 && numel(pOutIdx)==windowSamples+1 %due to rounding issues, in necessary, examine how this can be solved
                pOutIdx(1)=[];
            elseif pOutIdx(1)==0
                pOutIdx=pOutIdx+1;
            end
            %}
            
            for i=1:nCh
                data=zeros(size(pRecIdx),'int16')';currRec=1;
                for j=1:nWindows
                    if ~isempty(pSingleTrialTimeStamps{j})
                        fseek(obj.fidA(obj.n2sA(channels(i))),obj.headerSizeByte+(pSingleTrialTimeStamps{j}(1)-1)*obj.bytesPerRecCont+sum(obj.blkBytesCont(1:3)),'bof'); %(64+32+32+32)/8=20
                        data(:,(1:recordsPerTrial(j))+currRec-1)=fread(obj.fidA(obj.n2sA(channels(i))), [obj.dataSamplesPerRecord recordsPerTrial(j)], '1024*int16',obj.bytesPerRecCont - obj.blkBytesCont(4),'b');
                        currRec=currRec+recordsPerTrial(j);
                    else
                        disp('requested time stamp outside recording range!');
                    end
                end
                V_uV(pOutIdx+(i-1)*nWindows*windowSamples)=data(pRecIdx');
            end
            
            if obj.convertData2Double
                V_uV = permute(double(V_uV) * obj.MicrovoltsPerADAnalog(1) ,[3 2 1]);
            else
                V_uV = permute(V_uV,[3 2 1]);
            end
            
            if nargout==2
                t_ms=(1:windowSamples)*(1e3/obj.samplingFrequencyAnalog(1));
            end
        end
        
        function [T_ms,chNumber]=getTrigger(obj,startTime_ms,window_ms,name)
            %Extract triggers from file Neuralynx recording
            %Usage : [T_ms]=obj.getTrigger(startTime_ms,window_ms)
            %Input : name - which bit to extract for time stamps (out of 8,default = first bit, 1)
            %        startTime_ms - start time [ms].
            %        window_ms - the window duration [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %
            %Output: T_ms - trigger times [ms] - different triggers are arranged in a cell array
            
            %read Data from *.nev files (NeuroLynx) - files containing the data
            %fieldSelectionTrigger: (1): Timestamps, (2): Event IDs, (3): TTLs, (4): Extras, (5): Event Strings
            %Bits are coded as decimal integers
            
            timestamps = obj.segRead('timestamps')./obj.samplingFrequency(1)*1000-obj.globalStartTime_ms;
            %sampleNum = obj.segRead('sampleNum');
            eventType = obj.segRead('eventType');
            %nodeId = obj.segRead('nodeId');
            eventId = obj.segRead('eventId');
            ch = obj.segRead('data');
            activeCh=unique(ch);
            %if version >= 0.2, info.recNum = obj.segRead('recNum'); end
            
            if nargin<3
                pTime=true(size(eventType));
                startTime_ms=0;
            else
                pTime=timestamps>=startTime_ms & timestamps<(startTime_ms+window_ms);
            end
            timestamps=timestamps-startTime_ms;
            
            pTTL=eventType==3;
            if numel(activeCh)>0
            for i=1:numel(activeCh)
                T_ms{(2*activeCh(i)+1)}=timestamps(eventId==1 & ch==activeCh(i) & pTTL & pTime)';%ch1 is 0
                T_ms{(2*activeCh(i)+2)}=timestamps(eventId==0 & ch==activeCh(i) & pTTL & pTime)';%ch1 is 0
                chNumber((2*activeCh(i)+1))=activeCh(i);
                chNumber((2*activeCh(i)+2))=activeCh(i);
            end
            else
                disp('All event channels in the recording are empty');
                T_ms={};
                chNumber=[];
            end
        end



        function [camTrigs]=getCamerasTrigger(obj, trig_chanel, stop_len)
            %Extract triggers fron event channel, only for the times
            %between 2 triggers stops of above 1 sec.
            % It calls the regular getTiggers, and add another selection
            % for the relavnt triggers.
            % Usage : [camTrigs]=obj.getCamerasTrigger(startTime_ms,window_ms)
            %Input :
            % Output:camTrigs= triggers vectoer of times of on triggers, in ms from
            % the recording start time (0).
            if nargin < 3
                stop_len = 1000; % assuming that in Reptilearn the stop in the triggers are 2000 ms.
            end

            startTime_ms = 0;
            tTrig = getTrigger(obj,startTime_ms);
            trig = tTrig{trig_chanel}; % might change according to connections to the aquisition board.

            edges = find(diff(trig)>stop_len);
            if length(edges) == 2
                camTrigs = trig(edges(1)+1:edges(2)+1); % only the triggers between stops.
            elseif length(edges) >= 3
                disp('more than 2 stops in triggers, taking the triggers between the largest gaps.')

                maxD = find(max(diff(edges)));
                camTrigs = trig(edges(maxD)+1:edges(maxD+1)+1);

                %camTrigs = trig(edges(1)+1:edges(2)+1); % only the triggers between stops.
            elseif length(edges)==1
                disp('only 1 stop in triggers, check recording integrity.')
            end

        end

        function obj=closeOpenFiles(obj) %clear all open file handles
            if ~iscell(obj.recordingDir)
                nRecordings=1;
            else
                nRecordings=numel(obj.recordingDir);
            end
            if ~isempty(obj.fid)
                for j=1:nRecordings
                    fclose(obj.fid(j));
                end
            end
            if ~isempty(obj.fidA)
                for j=1:nRecordings
                    fclose(obj.fidA(j));
                end
            end
            if ~isempty(obj.fidEvnt)
                fclose(obj.fidEvnt);
            end
        end
        
        function obj=extractMetaData(obj)

            fprintf('\nExtracting meta data from: %s...',obj.recordingDir);

            channelFiles=dir([obj.recordingDir filesep '*.' obj.fileExtension]);
            channelFiles={channelFiles.name};
            nFiles=numel(channelFiles);

            %initialize struct with meta data for all channels
            fileData.fileSize=zeros(1,nFiles);
            fileData.samplingFrequency=zeros(1,nFiles);
            fileData.MicrovoltsPerAD=zeros(1,nFiles);
            fileData.startDate=cell(1,nFiles);
            fileData.blockLength=zeros(1,nFiles);
            fileData.dataDescriptionCont=cell(1,nFiles);
            fileData.fileHeaders=cell(1,nFiles);
            fileData.softwareVersion=zeros(1,nFiles);
            fileData.channelName=cell(1,nFiles);
            fileData.channelType=cell(1,nFiles);

            %go over files and extract meta data from each file
            for i=1:nFiles
                tmpFid=fopen([obj.recordingDir filesep channelFiles{i}],'r');
                if tmpFid==-1
                    error('The channel file identifier could not be initiated! Check that the experiment folder is correct and the files available for reading');
                end
                fseek(tmpFid,0,'eof');
                fileData.fileSize(i)=ftell(tmpFid);
                
                fseek(tmpFid,0,'bof');
                hdr = fread(tmpFid, obj.headerSizeByte, 'char*1');
                eval(char(hdr')); %convert to an header structure
                
                fileData.samplingFrequency(i)=header.sampleRate;
                fileData.MicrovoltsPerAD(i)=header.bitVolts;
                fileData.startDate{i}=header.date_created;
                %obj.bufferSize(i)=header.bufferSize(1);
                fileData.blockLength(i)=header.blockLength;
                fileData.dataDescriptionCont{i}=header.description;
                fileData.channelType{i}=header.channelType;
                fileData.fileHeaders{i} = header;
                if isfield(header, 'version')
                    fileData.softwareVersion(i) = header.version;
                    fileData.channelName{i}=header.channel;
                    fileData.channelType{i}=header.channelType;
                else
                    fileData.softwareVersion(i) = 0.0;
                end
                fclose(tmpFid);
            end

            %get data from xml file
            xmlData=[];
            if exist([obj.recordingDir filesep 'settings.xml'],'file')
                obj.openEphyXMLData = readstruct([obj.recordingDir filesep 'settings.xml']);
                
                if obj.openEphyXMLData.INFO.VERSION>"0.6.0"
                    obj.openEphyXMLStructureData = readstruct([obj.recordingDir filesep 'structure.openephys'],'FileType','xml');
                    xmlData.channelNamesAll={obj.openEphyXMLStructureData.RECORDING(1).STREAM.CHANNEL.nameAttribute};
                    xmlData.channelFiles=convertContainedStringsToChars({obj.openEphyXMLStructureData.RECORDING(1).STREAM.CHANNEL.filenameAttribute});
                    %settingBit2Volt={obj.openEphyXMLStructureData.RECORDING.STREAM.CHANNEL.bitVoltsAttribute};
                    
                    xmlData.chTypeName=cellfun(@(x) split(x,'_'),xmlData.channelNamesAll,'UniformOutput',0);
                    xmlData.chTypeName=cellfun(@(x) x(end),xmlData.chTypeName);
                    
                    xmlData.channelNames=cellfun(@(x) char(x),xmlData.channelNamesAll,'UniformOutput',0);
                    xmlData.channelNumbers=cellfun(@(x) str2double(regexp(x,'\d+','match')),xmlData.chTypeName);
                    xmlData.pAnalogCh=cellfun(@(x) x(1)=="A",cellfun(@(x) char(x),xmlData.chTypeName,'UniformOutput',0));
                    xmlData.pCh=cellfun(@(x) x(1)=="C",cellfun(@(x) char(x),xmlData.chTypeName,'UniformOutput',0));
                    obj.eventFileName=convertStringsToChars(obj.openEphyXMLStructureData.RECORDING(1).STREAM.EVENTS.filenameAttribute);
                else
                    %get channel information from settings files - including channel numbers - these may be in different order than the recorded files.
                    xmlData.processorNames={obj.openEphyXMLData.SIGNALCHAIN.PROCESSOR.pluginNameAttribute};
                    xmlData.pSources=find(cellfun(@(x) x=="Rhythm FPGA",xmlData.processorNames));
                    xmlData.chTypeName={obj.openEphyXMLData.SIGNALCHAIN.PROCESSOR(xmlData.pSources).CHANNEL_INFO.CHANNEL.nameAttribute};
                    xmlData.isRecorded={obj.openEphyXMLData.SIGNALCHAIN.PROCESSOR(xmlData.pSources).CHANNEL.SELECTIONSTATE};
                    xmlData.isRecorded=cellfun(@(x) x.paramAttribute,xmlData.isRecorded);
                    xmlData.isRecorded=logical(xmlData.isRecorded);

                    xmlData.channelNumbers=cellfun(@(x) str2double(regexp(x,'\d+','match')),xmlData.chTypeName);
                    xmlData.channelNames=xmlData.chTypeName(xmlData.isRecorded);
                    xmlData.channelNames=cellfun(@(x) char(x),xmlData.channelNames,'UniformOutput',false);%convert to adhere to the format of fileData
                    xmlData.pAnalogCh=cellfun(@(x) x(1)=="A",cellfun(@(x) char(x),xmlData.channelNames,'UniformOutput',0));
                    xmlData.pCh=cellfun(@(x) x(1)=="C",cellfun(@(x) char(x),xmlData.channelNames,'UniformOutput',0));
                end
            end

            %compare xlm data with channel data to verify consistency
            if ~isempty(xmlData)
                if numel(xmlData.chTypeName)==numel(fileData.channelName)
                    [C,ia,ib] = intersect(xmlData.channelNames,fileData.channelName);
                    %pAnalogCh=cellfun(@(x) x(1)=="A",cellfun(@(x) char(x),chTypeName(isRecorded),'UniformOutput',0));
                    %pCh=cellfun(@(x) x(1)=="C",cellfun(@(x) char(x),chTypeName(isRecorded),'UniformOutput',0));
                    if numel(C) == numel(fileData.channelName)
                        fprintf('xml data matches file data.\n')
                    end
                else
                    fprintf('There is mismatch between the settings xml and the data files!\nIf file data is correct, rename the xml file and run again. This will extract all information from the recording files\n');
                    %return;
                end
            else
                fprintf('settings.xml file not found! Extracting info from coninuous recording files.\n');
            end

            fileData.type=cellfun(@(x) regexp(x,'\S+','match'),fileData.channelName);
            fileData.channelNumbers=cellfun(@(x) str2double(regexp(x,'\d+$','match')),fileData.channelName);

            pAnalogCh=find(cellfun(@(x) x(1:3)=="ADC" || x(1:3)=="C1_",fileData.type));
            [~,p]=sort(fileData.channelNumbers(pAnalogCh));pAnalogCh=pAnalogCh(p);

            pCh=find(cellfun(@(x) x(1:2)=="CH",fileData.type));
            [~,p]=sort(fileData.channelNumbers(pCh));pCh=pCh(p);

            %{
            % for extracting events from file names in case this will be needed again
       

                channelFiles=dir([obj.recordingDir filesep '*.' obj.fileExtension]);
                channelFiles={channelFiles.name};
                dataFileStrings=cellfun(@(x) split(x(1:end-numel(obj.fileExtension)-1),'_'),channelFiles,'UniformOutput',0);
                nFields=numel(dataFileStrings{1});
                if nFields==2
                    recordNodes=cellfun(@(x) x{1},dataFileStrings,'UniformOutput',0);
                    channelNamesAll=cellfun(@(x) x{2},dataFileStrings,'UniformOutput',0);
                elseif nFields==3
                    recordNodes=cellfun(@(x) x{1},dataFileStrings,'UniformOutput',0);
                    recordNodeNames=cellfun(@(x) x{2},dataFileStrings,'UniformOutput',0);
                    channelNamesAll=cellfun(@(x) x{3},dataFileStrings,'UniformOutput',0);
                end
                channelNumbersAll=cellfun(@(x) str2double(regexp(x,'\d+','match')),channelNamesAll,'UniformOutput',1);
                channelPrefixAll=cellfun(@(x) regexp(x,'[A-Z]+','match'),channelNamesAll,'UniformOutput',1);
                pAnalogCh=cellfun(@(x) x(1:2)=="AU" | x(1:2)=="AD",cellfun(@(x) char(x),channelPrefixAll,'UniformOutput',0));%this may require updating for old versions of open ephys
                pCh=cellfun(@(x) x(1:2)=="CH",cellfun(@(x) char(x),channelPrefixAll,'UniformOutput',0));%this may require updating for old versions of open ephys
            
                           %get channel information from files in folder
                    obj.eventFileName='all_channels.events';
                    channelFiles=dir([obj.recordingDir filesep '*.' obj.fileExtension]);
                    channelFiles={channelFiles.name};
                    dataFileStrings=cellfun(@(x) split(x(1:end-numel(obj.fileExtension)-1),'_'),channelFiles,'UniformOutput',0);
                    nFields=numel(dataFileStrings{1});
                    if nFields==2
                        recordNodesFromFile=cellfun(@(x) x{1},dataFileStrings,'UniformOutput',0);
                        channelNamesFromFile=cellfun(@(x) x{2},dataFileStrings,'UniformOutput',0);
                    elseif nFields==3
                        recordNodesFromFile=cellfun(@(x) x{1},dataFileStrings,'UniformOutput',0);
                        recordNodeNamesFromFile=cellfun(@(x) x{2},dataFileStrings,'UniformOutput',0);
                        channelNamesFromFile=cellfun(@(x) x{3},dataFileStrings,'UniformOutput',0);
                    end
                    try
                        channelNumbersFromFile=cellfun(@(x) str2double(regexp(x,'\d+','match')),channelNamesFromFile,'UniformOutput',1);
                        chTypeNameFromFile=cellfun(@(x) regexp(x,'[A-Z]+','match'),channelNamesFromFile,'UniformOutput',1);
                        pAnalogChFromFile=cellfun(@(x) x(1)=="A",cellfun(@(x) char(x),chTypeNameFromFile,'UniformOutput',0));
                        pChFromFile=cellfun(@(x) x(1)=="C",cellfun(@(x) char(x),chTypeNameFromFile,'UniformOutput',0));
                        channelInfoInFileName=true;
                    catch
                        channelInfoInFileName=false;
                    end

            %}

            obj.channelFilesAnalog=channelFiles(pAnalogCh);
            obj.channelFiles=channelFiles(pCh);

            obj.channelNumbers=fileData.channelNumbers(pCh);
            if any(obj.channelNumbers==0)
                obj.channelNumbers=obj.channelNumbers+1;
            end
            obj.analogChannelNumbers=fileData.channelNumbers(pAnalogCh);
            if any(obj.analogChannelNumbers==0)
                obj.analogChannelNumbers=obj.analogChannelNumbers+1;
            end            
            obj.channelNames=fileData.channelName(pCh);
            obj.analogChannelNames=fileData.channelName(pAnalogCh);
            
            [obj.channelNumbers,pTmp]=sort(obj.channelNumbers);
            obj.channelFiles=obj.channelFiles(pTmp);
            obj.channelNames=obj.channelNames(pTmp);
            obj.n2s(obj.channelNumbers)=1:numel(obj.channelNumbers);
            
            nAnalogChannels=numel(obj.analogChannelNumbers);
            if numel(unique(obj.analogChannelNumbers))~=nAnalogChannels
                obj.analogChannelNumbers=1:nAnalogChannels;
                pTmp=1:nAnalogChannels;
                fprintf('Analog channel numbers contain duplicates!!! Reordering numbers serially.\n');
            else
                [obj.analogChannelNumbers,pTmp]=sort(obj.analogChannelNumbers);
            end
            obj.channelFilesAnalog=obj.channelFilesAnalog(pTmp);
            obj.analogChannelNames=obj.analogChannelNames(pTmp);
            obj.n2sA(obj.analogChannelNumbers)=1:nAnalogChannels;

            %fill all properties with the relevant keywords
            obj.MicrovoltsPerAD=fileData.MicrovoltsPerAD(pCh);
            obj.MicrovoltsPerADAnalog=fileData.MicrovoltsPerAD(pAnalogCh)*1e6; %in open ephys the AD value for AUX channels is given in volts.
            obj.ZeroADValue=zeros(size(obj.MicrovoltsPerAD));
            obj.ZeroADValueAnalog=zeros(size(obj.MicrovoltsPerADAnalog));
            obj.samplingFrequency=fileData.samplingFrequency(pCh);
            obj.samplingFrequencyAnalog=fileData.samplingFrequency(pAnalogCh);
            obj.fileSize=fileData.fileSize(pCh);
            obj.fileSizeAnalog=fileData.fileSize(pAnalogCh);
            obj.softwareVersion=fileData.softwareVersion(pCh);
            obj.softwareVersionAnalog=fileData.softwareVersion(pAnalogCh);
            obj.blockLength=fileData.blockLength(pCh);
            obj.blockLengthAnalog=fileData.blockLength(pAnalogCh);
            obj.startDate=fileData.startDate(pCh);
            obj.startDateA=fileData.startDate(pAnalogCh);
            obj.fileHeaders=fileData.fileHeaders(pCh);
            obj.fileHeadersAnalog=fileData.fileHeaders(pAnalogCh);
            obj.dataDescriptionCont=fileData.dataDescriptionCont(pCh);
            obj.dataDescriptionContAnalog=fileData.dataDescriptionCont(pAnalogCh);

            obj=obj.getFileIdentifiers;

            %prepare data structures for continuous - assumes that all channels have the same size structure and time stamps, else run separetly on every file
            bStr = {'ts' 'nsamples' 'recNum' 'data' 'recordMarker'};
            bTypes = {'int64' 'uint16' 'uint16' 'int16' 'uint8'};
            bRepeat = {1 1 1 obj.dataSamplesPerRecord 10};
            obj.blkCont = struct('Repeat',bRepeat,'Types', bTypes,'Str',bStr);
            if obj.softwareVersion < 0.2, obj.blkCont(3) = []; end
            if obj.softwareVersion < 0.1, obj.blkCont(1).Types = 'uint64'; obj.blkCont(2).Types = 'int16'; end
            obj.blkBytesCont = str2double(regexp({obj.blkCont.Types},'\d{1,2}$','match', 'once')) ./8 .* cell2mat({obj.blkCont.Repeat});
            obj.bytesPerRecCont=sum(obj.blkBytesCont);
            
            fprintf('\nExtracting time stamp information...');
            if numel(obj.channelNumbers)==0
                disp('0 electrode channels found in recording!');
                obj.samplingFrequency(1)=obj.samplingFrequencyAnalog(1);
                obj.nRecordsCont = floor((obj.fileSizeAnalog(1) - obj.headerSizeByte)/obj.bytesPerRecCont);
                fid4Tests=obj.fidA(1);
            else
                obj.nRecordsCont = floor((obj.fileSize(1) - obj.headerSizeByte)/obj.bytesPerRecCont);
                fid4Tests=obj.fid(1);
            end
            obj.recordLength = obj.dataSamplesPerRecord/obj.samplingFrequency(1)*1000; %must appear after extracting sampling frequency.
            obj.sample_ms=1e3/obj.samplingFrequency(1);

            fseek(fid4Tests, obj.headerSizeByte, 'bof');
            obj.allTimeStamps = fread(fid4Tests, obj.nRecordsCont*obj.blkCont(1).Repeat, sprintf('%d*%s', obj.blkCont(1).Repeat,obj.blkCont(1).Types), obj.bytesPerRecCont - obj.blkBytesCont(1), 'l')/obj.samplingFrequency(1)*1000;
            obj.globalStartTime_ms=obj.allTimeStamps(1);
            obj.allTimeStamps = obj.allTimeStamps-obj.globalStartTime_ms; %time stamp vector must be a column vector!!!
            if obj.allTimeStamps(end)<=obj.allTimeStamps(end-1)
                obj.allTimeStamps(end)=[];
                fprintf('Last record has a non-valid stims stamp!!! Removing last data record.\n');
            end
            obj.recordingDuration_ms=obj.allTimeStamps(end);
            if any(diff(obj.allTimeStamps)>obj.dataSamplesPerRecord), fprintf(2,'\nError!!! Some blocks are missing in recording!!!\n'),end
            
            %check that all records have 1024 samples (data integrity check)
            fprintf('\nChecking integrity of all records in ch1...');
            fseek(fid4Tests, obj.headerSizeByte+obj.blkBytesCont(1), 'bof');
            sampleNumbers = fread(fid4Tests, obj.nRecordsCont*obj.blkCont(2).Repeat, sprintf('%d*%s', obj.blkCont(2).Repeat,obj.blkCont(2).Types), obj.bytesPerRecCont - obj.blkBytesCont(2), 'l');
            if ~all(sampleNumbers == obj.dataSamplesPerRecord) && obj.version >= 0.1, error('Found currupted records!!! Please check manually'); end
            
            %prepare data structures for events
            bStr = {'timestamps' 'sampleNum' 'eventType' 'nodeId' 'eventId' 'data' 'recNum'};
            bTypes = {'int64' 'uint16' 'uint8' 'uint8' 'uint8' 'uint8' 'uint16'};
            bRepeat = {1 1 1 1 1 1 1};
            obj.blkEvnt = struct('Repeat',bRepeat,'Types', bTypes,'Str',bStr);
            if obj.softwareVersion < 0.2, obj.blkEvnt(7) = [];  end
            if obj.softwareVersion < 0.1, obj.blkEvnt(1).Types = 'uint64'; end
            
            fseek(obj.fidEvnt,0,'eof');
            obj.evntFileSize=ftell(obj.fidEvnt);
            
            obj.blkBytesEvnt = str2double(regexp({obj.blkEvnt.Types},'\d{1,2}$','match', 'once')) ./8 .* cell2mat({obj.blkEvnt.Repeat});
            obj.bytesPerRecEvnt=sum(obj.blkBytesEvnt);
            obj.nRecordsEvnt = floor((obj.evntFileSize - obj.headerSizeByte)/obj.bytesPerRecEvnt);
            
            %{
            %prepare data structures for spikes
            num_channels = info.header.num_channels;
            num_samples = 40;
            bStr = {'eventType' 'timestamps' 'timestamps_software' 'source' 'nChannels' 'nSamples' 'sortedId' 'electrodeID' 'channel' 'color' 'pcProj' 'samplingFrequencyHz' 'data' 'gain' 'threshold' 'recordingNumber'};
            bTypes = {'uint8' 'int64' 'int64' 'uint16' 'uint16' 'uint16' 'uint16' 'uint16' 'uint16' 'uint8' 'float32' 'uint16' 'uint16' 'float32' 'uint16' 'uint16'};
            bRepeat = {1 1 1 1 1 1 1 1 1 3 2 1 num_channels*num_samples num_channels num_channels 1};
            obj.blkSpk = struct('Repeat',bRepeat,'Types', bTypes,'Str',bStr);
            if obj.softwareVersion < 0.4,  obj.blkSpk(7:12) = []; obj.blkSpk(8).Types = 'uint16'; end
            if obj.softwareVersion == 0.3, obj.blkSpk = [obj.blkSpk(1), struct('Repeat',1,'Types','uint32','Str','ts'), obj.blkSpk(2:end)]; end
            if obj.softwareVersion < 0.3, obj.blkSpk(2) = []; end
            if obj.softwareVersion < 0.2, obj.blkSpk(9) = []; end
            if obj.softwareVersion < 0.1, obj.blkSpk(2).Types = 'uint64'; end
            
            %}
            
            fprintf('\nsaving meta data\n');
            obj.saveMetaData;
        end
        
        function [obj]=getFileIdentifiers(obj)
            %get file indentifiers for all channel files
            for i=1:numel(obj.channelFiles)
                obj.fid(i)=fopen([obj.recordingDir filesep obj.channelFiles{i}],'r');
            end
            if any(obj.fid==-1)
                error('The channel file identifier could not be created. Check that the channel file names are correct and available for reading');
            end
            for i=1:numel(obj.channelFilesAnalog)
                obj.fidA(i)=fopen([obj.recordingDir filesep obj.channelFilesAnalog{i}],'r');
            end
            if any(obj.fidA==-1)
                error('The analog file identifier could not be created. Check that the analog file names are correct and available for reading');
            end
            if ~isempty(obj.eventFileName)
                obj.fidEvnt=fopen([obj.recordingDir filesep obj.eventFileName],'r');
            end
            if any(obj.fidEvnt==-1)
                error('The event file identifier could not be created. Check that the event file names are correct and available for reading');
            end
        end
        
        function seg = segRead(obj,segName, mf)
            if nargin == 2, mf = 'l'; end
            segNum = find(strcmp({obj.blkEvnt.Str},segName));
            fseek(obj.fidEvnt, sum(obj.blkBytesEvnt(1:segNum-1))+obj.headerSizeByte, 'bof');
            seg = fread(obj.fidEvnt, obj.nRecordsEvnt*obj.blkEvnt(segNum).Repeat, sprintf('%d*%s', obj.blkEvnt(segNum).Repeat,obj.blkEvnt(segNum).Types), obj.bytesPerRecEvnt - obj.blkBytesEvnt(segNum), mf);
        end
        
        function obj=getRecordingFiles(obj,recordingFile,fileExtension)
            %Get directory with data files
            %Usage: obj = getRecordingFiles(obj,recordingFile,fileExtension)
            %if no recording file is entered lauches a GUI
            %if no file extension entered, a directory is chosen rather than a specific files (for example for neuralynx recordings)
            obj.folderMode=1;
            obj.multifileMode=false;
            
            %If no files were entered open GUI for choosing a file or a directory else get the files entered
            if ~isempty(recordingFile) %if directory with data was not entered open get directory GUI
                
                if obj.folderMode
                    if recordingFile(end)~=filesep
                        recordingFile=[recordingFile filesep];
                    end
                end
                [pathstr, name, ext] = fileparts(recordingFile);
                
                
                if ~exist([pathstr filesep obj.dataFileNames],'dir')
                    disp(['Searching for recording file: ' [pathstr filesep obj.dataFileNames]]);
                    error('Object was not constructed since no valid recording file name was chosen');
                end

                obj.dataFileNames=recordingFile;
                
            else %if directory with data was not entered open get directory GUI
                [obj.recordingDir]= uigetdir(obj.defaultLocalDir,'Choose the data folder');
                [pathstr, name] = fileparts(obj.recordingDir);
                
            end
            obj.recordingDir=[pathstr filesep name];
            obj.recordingName = name;
            obj.metaDataFile=[obj.recordingDir filesep obj.recordingName 'OE_metaData.mat'];
            %obj.dataFileNames=dir([pathstr filesep name filesep '*.' obj.fileExtension]);
            
            if obj.recordingDir(end)=='/'
                obj.recordingDir=obj.recordingDir(1:end-1);
            end
        end

    end
    
    methods (Hidden)
        %class constructor
        function obj = OERecording(recordingFile)
            %Usage: obj = NeuraLynxRecording(recordingFile)
            obj.datatype='int16';
            obj.folderMode=1;
            %get data files
            if nargin==0
                recordingFile=[];
            elseif nargin>1
                disp('Object was not constructed since too many parameters were given at construction');
                return;
            end
            if iscell(recordingFile)
                recordingFile=recordingFile{1};
            end
            obj=obj.getRecordingFiles(recordingFile);
                        
            if exist(obj.metaDataFile,'file') && ~obj.overwriteMetaData
                obj=loadMetaData(obj);
                obj=obj.getFileIdentifiers;
            else
                obj=extractMetaData(obj);%file identifiers are aquired inside function
            end
            
            %load layout
            obj=obj.loadChLayout;
            if isempty(obj.chLayoutNumbers)
                disp('No layout file in directory, to add layout please add a *.chMap file to the data directory with the layout name');
            end
            
        end
        
    end
    
end

%{
function [data, timestamps, info] = load_open_ephys_data_faster(filename)
%
% [data, timestamps, info] = load_open_ephys_data(filename)
%
%   Loads continuous, event, or spike data files into Matlab.
%
%   Inputs:
%
%     filename: path to file
%
%
%   Outputs:
%
%     data: either an array continuous samples (in microvolts),
%           a matrix of spike waveforms (in microvolts),
%           or an array of event channels (integers)
%
%     timestamps: in seconds
%
%     info: structure with header and other information
%
%
%
%   DISCLAIMER:
%
%   Both the Open Ephys data format and this m-file are works in progress.
%   There's no guarantee that they will preserve the integrity of your
%   data. They will both be updated rather frequently, so try to use the
%   most recent version of this file, if possible.
%
%

%
%     ------------------------------------------------------------------
%
%     Copyright (C) 2014 Open Ephys
%
%     ------------------------------------------------------------------
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     <http://www.gnu.org/licenses/>.
%

[~,~,filetype] = fileparts(filename);
if ~any(strcmp(filetype,{'.events','.continuous','.spikes'}))
    error('File extension not recognized. Please use a ''.continuous'', ''.spikes'', or ''.events'' file.');
end

fid = fopen(filename);
fseek(fid,0,'eof');
filesize = ftell(fid);

NUM_HEADER_BYTES = 1024;
fseek(fid,0,'bof');
hdr = fread(fid, NUM_HEADER_BYTES, 'char*1');
info = getHeader(hdr);
if isfield(info.header, 'version')
    version = info.header.version;
else
    version = 0.0;
end

switch filetype
    case '.events'
        bStr = {'timestamps' 'sampleNum' 'eventType' 'nodeId' 'eventId' 'data' 'recNum'};
        bTypes = {'int64' 'uint16' 'uint8' 'uint8' 'uint8' 'uint8' 'uint16'};
        bRepeat = {1 1 1 1 1 1 1};
        dblock = struct('Repeat',bRepeat,'Types', bTypes,'Str',bStr);
        if version < 0.2, dblock(7) = [];  end
        if version < 0.1, dblock(1).Types = 'uint64'; end
    case '.continuous'
        SAMPLES_PER_RECORD = 1024;
        bStr = {'ts' 'nsamples' 'recNum' 'data' 'recordMarker'};
        bTypes = {'int64' 'uint16' 'uint16' 'int16' 'uint8'};
        bRepeat = {1 1 1 SAMPLES_PER_RECORD 10};
        dblock = struct('Repeat',bRepeat,'Types', bTypes,'Str',bStr);
        if version < 0.2, dblock(3) = []; end
        if version < 0.1, dblock(1).Types = 'uint64'; dblock(2).Types = 'int16'; end
    case '.spikes'
        num_channels = info.header.num_channels;
        num_samples = 40;
        bStr = {'eventType' 'timestamps' 'timestamps_software' 'source' 'nChannels' 'nSamples' 'sortedId' 'electrodeID' 'channel' 'color' 'pcProj' 'samplingFrequencyHz' 'data' 'gain' 'threshold' 'recordingNumber'};
        bTypes = {'uint8' 'int64' 'int64' 'uint16' 'uint16' 'uint16' 'uint16' 'uint16' 'uint16' 'uint8' 'float32' 'uint16' 'uint16' 'float32' 'uint16' 'uint16'};
        bRepeat = {1 1 1 1 1 1 1 1 1 3 2 1 num_channels*num_samples num_channels num_channels 1};
        dblock = struct('Repeat',bRepeat,'Types', bTypes,'Str',bStr);
        if version < 0.4,  dblock(7:12) = []; dblock(8).Types = 'uint16'; end
        if version == 0.3, dblock = [dblock(1), struct('Repeat',1,'Types','uint32','Str','ts'), dblock(2:end)]; end
        if version < 0.3, dblock(2) = []; end
        if version < 0.2, dblock(9) = []; end
        if version < 0.1, dblock(2).Types = 'uint64'; end
end
blockBytes = str2double(regexp({dblock.Types},'\d{1,2}$','match', 'once')) ./8 .* cell2mat({dblock.Repeat});
numIdx = floor((filesize - NUM_HEADER_BYTES)/sum(blockBytes));

switch filetype
    case '.events'
        timestamps = segRead('timestamps')./info.header.sampleRate;
        info.sampleNum = segRead('sampleNum');
        info.eventType = segRead('eventType');
        info.nodeId = segRead('nodeId');
        info.eventId = segRead('eventId');
        data = segRead('data');
        if version >= 0.2, info.recNum = segRead('recNum'); end
    case '.continuous'
        info.ts = segRead('ts');
        info.nsamples = segRead('nsamples');
        if ~all(info.nsamples == SAMPLES_PER_RECORD)&& version >= 0.1, error('Found corrupted record'); end
        if version >= 0.2, info.recNum = segRead('recNum'); end
        data = segRead('data', 'b').*info.header.bitVolts; % read in data
        timestamps = nan(size(data));
        current_sample = 0;
        for record = 1:length(info.ts)
            timestamps(current_sample+1:current_sample+info.nsamples(record)) = info.ts(record):info.ts(record)+info.nsamples(record)-1;
            current_sample = current_sample + info.nsamples(record);
        end
    case '.spikes'
        timestamps = segRead('timestamps')./info.header.sampleRate;
        info.source = segRead('source');
        info.samplenum = segRead('nSamples');
        info.gain = permute(reshape(segRead('gain'), num_channels, numIdx), [2 1]);
        info.thresh = permute(reshape(segRead('threshold'), num_channels, numIdx), [2 1]);
        if version >= 0.4, info.sortedId = segRead('sortedId'); end
        if version >= 0.2, info.recNum = segRead('recordingNumber'); end
        data = permute(reshape(segRead('data'), num_samples, num_channels, numIdx), [3 1 2]);
        data = (data-32768)./ permute(repmat(info.gain/1000,[1 1 num_samples]), [1 3 2]);
end
fclose(fid);

function seg = segRead(segName, mf)
    if nargin == 1, mf = 'l'; end
    segNum = find(strcmp({dblock.Str},segName));
    fseek(fid, sum(blockBytes(1:segNum-1))+NUM_HEADER_BYTES, 'bof');
    seg = fread(fid, numIdx*dblock(segNum).Repeat, sprintf('%d*%s', dblock(segNum).Repeat,dblock(segNum).Types), sum(blockBytes) - blockBytes(segNum), mf);
end

end
function info = getHeader(hdr)
eval(char(hdr'));
info.header = header;
end

%}