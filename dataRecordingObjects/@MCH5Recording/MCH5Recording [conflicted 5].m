classdef MCH5Recording < dataRecording
  %MCH5Recording Reads a h5 recording file exported by MCS DataManager.
  % This class is used to read h5 containing raw data and
  % trigger information (currently from digital in data) from a 
  % MultiChannelSystems h5 file exported by MultiChannel DataManager.
  % .h5 file specification can be found in the link:
  % https://mcspydatatools.readthedocs.io/en/latest/mcs_hdf5_protocol.html

  % AUTHOR: stefano.masneri@brain.mpg.de & Mark Shein-Idelson
  % DATE: 18.04.2017
  % 
  % TODO:
  %   - See how multiple recoding in a single h5 file works
  %   - Set up protocol for using triggers (i.e. digital events) instead of
  %   using get_triggers_from_digin. Maybe use check flag as "digital_events=1" for when
  %   it is set up vs "digital_events=0" for when we have to calc from dig in
  
  properties
    numRecordings;   % number of recordings in a single .h5 file. Currently only supports 1
    timestamps;      % timestamps information for each channel
%     triggerFilename; % name of the *.kwe file containing trigger info
    digital_events=0;
%     bitDepth;        % number of bits used to store data
    sample_ms;
    fullFilename;    % path + name
    recNameHD5;  % names of all recordings
    dataLength;      % total samples in the data
    totalChannels    % total channels in the data
    info;            % information on recording file
    analogInfo       % information on the analog stream of first recording
    lengthInfo;      % information on recording file length
%     globalStartTime  % start time within the session
    streamPaths      % paths to streams
    streamsSubTypes  % type pf streams (electrode, aux, digital)
    electrodeStreamNum %the stream number in h5 file containing raw electrode data (usually 0)
    auxStreamNum
    digitalStreamNum
    pathToRawDataStreamGroup
    pathToDigitalDataStreamGroup
    pathToAuxDataStreamGroup
    electrodeInfoChannel
%     digitalInfoChannel
%     auxInfoChannel
    unit %Physical unit of the measured sensor value
    exponent %1xn
    
  end
  
   properties (Constant, Hidden)
        defaultRawDataStreamName='Electrode';
        defaultFilteredDataStreamName='Filtere';
        defaultAnalogDataStreamName='Auxiliary';
        defaultTriggerStreamName='Trigger'; %CHANGE THIS TO EVENTS
        defaultDigitalDataStreamName='Digital';
%         fileExtension='h5';
        pathToAllRecordings='/Data/';
        pathToRecording='/Data/Recording_0/';
        pathToAnalogStream = '/Data/Recording_0/AnalogStream/'; %This is in case there is only one recording in HDF5 file. If there are more - code should be revised
   end
    
  properties (Constant = true)
    
    
%     pathToTriggerData = '/event_types/TTL/events/time_samples'; % where the triggers are stored in the .kwe file
%     pathToTriggerOnOff = '/event_types/TTL/events/user_data/eventID';
%     pathToTriggerChannel = '/event_types/TTL/events/user_data/event_channels';
        
    %must define them, bc abstract in base class
    defaultLocalDir='E:\Yuval\DataAnalysis'; %Default directory from which search starts
    signalBits = 24 %the quantization of the sampling card
    numberOfCharFromEndToBaseName=7;
  end
  
  methods
    
    function [V_uV ,t_ms]=getData(obj, channels, startTime_ms, window_ms)
      %GETDATA Extract MCS h5 recording data from file to memory
      %Usage: [V_uV,t_ms] = obj.getData(channels,startTime_ms,window_ms);
      %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
      %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
      %        window_ms - a scalar [1x1] with the window duration [ms].
      %Output: V_uv - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
      %        t_ms - A time vector relative to recording start (t=0 at start)
      
      windowSamples = double(round(double(window_ms) / obj.sample_ms(1)));
      nWindows = numel(startTime_ms);
      startTime_ms = round(startTime_ms/obj.sample_ms(1))*obj.sample_ms(1);
      startElement = double(round(startTime_ms/obj.sample_ms(1)));
      if startElement(1) == 0 || Inf == startElement(1)
        startElement(1) = 1;
      end
      %window_ms = windowSamples * obj.sample_ms(1);
      
      if isempty(channels) %if no channels are entered, get all channels
        channels=obj.channelNumbers;
      end
      channels=obj.n2s(channels); %check with mark
      nCh = numel(channels);
      
      V_uV = zeros(nCh, nWindows, windowSamples, obj.datatype); %initialize waveform matrix
      %{
      	YUVAL: I have transposed all channel/samples indices because for
      	some reason the data matrix is fliped here (i.e. when loading the
      	whole data into a matrix we get samplesXchannels matrix instead of
      	a channelsXsamples matrix.
        So check to see if this really works and/or if is consistent with
        larger datasets
        I have also transposed the h5read returned matrix
      %}
      % Speed up if all channels are consecutives
      if all(diff(channels)==1)  
        for m = 1:numel(startElement)
          if startElement(m) <= -windowSamples
            %do nothing, return all zeros
          elseif startElement(m) < 1
            V_uV(:, m, -startElement(m)+1 : end) = h5read(obj.fullFilename, [obj.pathToRawDataStreamGroup '/ChannelData'], ...
              [channels(1) 1]', [length(channels) windowSamples + startElement(m)]')';
          elseif startElement(m) >= obj.dataLength
            startElement(m) = obj.dataLength - windowSamples;
            V_uV(:, m, :) = h5read(obj.fullFilename, [obj.pathToRawDataStreamGroup '/ChannelData'], ...
              [channels(1) startElement(m)]', [length(channels) windowSamples]')';
          elseif startElement(m) + windowSamples > obj.dataLength
            V_uV(:, m, 1:obj.dataLength-startElement(m)) = h5read(obj.fullFilename, [obj.pathToRawDataStreamGroup '/ChannelData'], ...
              [channels(1) startElement(m)]', [length(channels) obj.dataLength - startElement(m)]')';
          else
            V_uV(:, m, :) = h5read(obj.fullFilename, [obj.pathToRawDataStreamGroup '/ChannelData'], ...
              [channels(1) startElement(m)]', [length(channels) windowSamples]')';
          end
        end
      else
        for k = 1:length(channels)
          for m = 1:numel(startElement)

            if startElement(m) <= -windowSamples
              %do nothing, return all zeros
            elseif startElement(m) < 1
              V_uV(k, m, -startElement(m)+1 : end) = h5read(obj.fullFilename, [obj.pathToRawDataStreamGroup '/ChannelData'], ...
                [channels(k) 1]', [1 windowSamples + startElement(m)]')';
            elseif startElement(m) >= obj.dataLength
              startElement(m) = obj.dataLength - windowSamples;
              V_uV(k, m, :) = h5read(obj.fullFilename, [obj.pathToRawDataStreamGroup '/ChannelData'], ...
                [channels(k) startElement(m)]', [1 windowSamples]')';
            elseif startElement(m) + windowSamples > obj.dataLength
              V_uV(k, m, 1:obj.dataLength-startElement(m)) = h5read(obj.fullFilename, [obj.pathToRawDataStreamGroup '/ChannelData'], ...
                [channels(k) startElement(m)]', [1 obj.dataLength - startElement(m)]')';
            else
              V_uV(k, m, :) = h5read(obj.fullFilename, [obj.pathToRawDataStreamGroup '/ChannelData'], ...
                [channels(k) startElement(m)]', [1 windowSamples]')';
            end
          end
        end
      end
      
      
      
      if obj.convertData2Double
          V_uV=double(V_uV);
          for k = 1:size(V_uV, 1)
              V_uV(k, :, :) = (V_uV(k, :, :)-obj.ZeroADValue(k)) * obj.MicrovoltsPerAD(k)*10^(obj.Exponent(k)+6); %exponent brings value in V, we want uV
          end
      end
      
      if nargout==2
        t_ms=(1:windowSamples)*(1e3/obj.samplingFrequency(1));
      end

    end
    
    function [T_ms, chNumber] = getTrigger(obj, ~, ~, ~)
      %GETTRIGGER Extract triggers from file Kwik recording
      % Read Data from *.kwe files (NeuroLynx) - files containing the data
      %Usage : [T_ms] = obj.getTrigger(name,startTime_ms,endTime_ms)
      %Input : name - which bit to extract for time stamps (out of 8,default = first bit, 1)
      %        startTime_ms - start time [ms].
      %        window_ms - the window duration [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
      %
      %Output: T_ms - trigger times [ms] - different triggers are arranged in a cell array
      
      allTriggers = h5read(obj.triggerFilename, obj.pathToTriggerData);
      onOff =  h5read(obj.triggerFilename, obj.pathToTriggerOnOff);
      %T_ms{1} = allTriggers(onOff == 1);
      %T_ms{2} = allTriggers(onOff == 0);
      chNumber = h5read(obj.triggerFilename, obj.pathToTriggerChannel);
      nTrig=unique(chNumber);
      T_ms = cell(1, numel(nTrig)*2);
      for i=1:numel(nTrig)
          pCh=chNumber==nTrig(i);
          T_ms{i*2-1}=double(allTriggers(onOff == 1 & pCh)-obj.globalStartTime)/(obj.samplingFrequency(1)/1000);
          T_ms{i*2}=double(allTriggers(onOff == 0 & pCh)-obj.globalStartTime)/(obj.samplingFrequency(1)/1000);
      end
    end
  end
  
  methods (Hidden = true)
    
    %class constructor
    function obj = MCH5Recording(recordingFile)
      %get data files
      if nargin == 0
        recordingFile=[];
      elseif nargin>1
        disp('MCH5Recording: Object was not constructed since too many parameters were given at construction');
        return;
      end
%       obj.datatype='int32';

      obj = obj.getRecordingFiles(recordingFile, 'h5');
      
      % Find the .h5 file
      %[~, name, ~] = fileparts(obj.recordingName);
      obj.fullFilename = fullfile(obj.recordingDir, obj.recordingName);
      filePrefix=strsplit(obj.recordingName,'_');
      
%       %%%YUVAL: CHANGE TO GET RELEVANT TRIGGER PROPERTIES FROM HDF5 %%%%%
%       triggerFile = dir([obj.recordingDir filesep filePrefix{1} '*.kwe']);
%       if isempty(triggerFile)
%           triggerFile = dir([obj.recordingDir filesep '*.kwe']);
%           disp(['Trigger file with prefix ' filePrefix{1} ' not found!!!!, looking for other .kwe files in the same folder']);
%       end
%       if isempty(triggerFile)
%         error('KwikRecording: Cannot file .kwe file')
%       elseif length(triggerFile) > 1
%         warning('KwikRecording: Multiple .kwe file found! using the first one')
%       end
%       obj.triggerFilename = fullfile(obj.recordingDir, triggerFile.name);
      
      if exist([obj.recordingDir filesep 'metaData.mat'],'file') && ~obj.overwriteMetaData
          obj = loadMetaData(obj); %needs recNameHD5
      else
          obj = extractMetaData(obj);
      end
      
      obj.numRecordings = length(obj.info.Groups);
      if obj.numRecordings > 1
          warning('KwikRecording: file contains multiple recordings.')
      end
      
      obj.recNameHD5 = cell(1, obj.numRecordings);
      for k = 1:obj.numRecordings
          obj.recNameHD5{k} = obj.info.Groups(k).Name;
      end

     obj.dataLength = obj.lengthInfo.Dataspace.Size(1);
     obj.totalChannels=obj.lengthInfo.Dataspace.Size(2);
     
      %layout
      obj=obj.loadChLayout;
      if isempty(obj.chLayoutNumbers)
          if obj.totalChannels<=32
              obj.layoutName='layout_300_6x6_FlexMEA';
              obj.electrodePitch=300;
          elseif obj.totalChannels<=60
              obj.layoutName='layout_200_8x8.mat';
              obj.electrodePitch=200;
          elseif obj.totalChannels<=120
              obj.layoutName='layout_100_12x12.mat';
              obj.electrodePitch=100;
          elseif obj.totalChannels<=252
              obj.layoutName='layout_100_16x16.mat';
              obj.electrodePitch=100;
          end
          load(obj.layoutName);
          obj.chLayoutNumbers=En;
          obj.chLayoutNames=Ena;
          obj.chLayoutPositions=Enp;
      end
      
    end
    
    function delete(obj) %do nothing
    end
    
  end
  
  methods
    
    function obj = extractMetaData(obj)
        %a lot of these data come as a vector specifing a value for each \
        %channel. assuming they are the same for all, just the first one 
        %is taken (i.e. time tick, ADZero,RawDataType,ConcersionFactor,
        %HighPassFilterType etc
        if ~strcmp(obj.fullFilename(end-1:end),'h5')
            obj.fullFilename=[obj.fullFilename '.h5'];
        end
        obj.info=h5info(obj.fullFilename, obj.pathToAllRecordings);
        obj.analogInfo = h5info(obj.fullFilename, obj.pathToAnalogStream);
        obj.streamPaths{1}=obj.analogInfo.Groups(1).Name; %should be stream 0 path
        obj.streamPaths{2}=obj.analogInfo.Groups(2).Name; %should be stream 1 path
        obj.streamPaths{3}=obj.analogInfo.Groups(3).Name; %should be stream 2 path
        
        %get streams' numbers and paths (for raw electrode, digital and aux
        %streams. Usually it's 0,1,2 but just to make sure
        obj.streamsSubTypes={h5readatt(obj.fullFilename,obj.streamPaths{1},'DataSubType'),h5readatt(obj.fullFilename,obj.streamPaths{2},'DataSubType'),h5readatt(obj.fullFilename,obj.streamPaths{3},'DataSubType')};
        obj.electrodeStreamNum=find(ismember(obj.streamsSubTypes,obj.defaultRawDataStreamName))-1; %'Electrode';'Auxiliary';'Digital';
        obj.auxStreamNum=find(ismember(obj.streamsSubTypes,obj.defaultAnalogDataStreamName))-1;
        obj.digitalStreamNum=find(ismember(obj.streamsSubTypes,obj.defaultDigitalDataStreamName))-1;
       
        obj.pathToRawDataStreamGroup=obj.streamPaths{obj.electrodeStreamNum+1}; %these are without '/' ending 
        obj.pathToDigitalDataStreamGroup=obj.streamPaths{obj.digitalStreamNum+1};
        obj.pathToAuxDataStreamGroup=obj.streamPaths{obj.auxStreamNum+1};
        
        obj.lengthInfo = h5info(obj.fullFilename, [obj.pathToRawDataStreamGroup '/ChannelData']);
        
        obj.electrodeInfoChannel=h5read(obj.fullFilename, [obj.pathToRawDataStreamGroup '/InfoChannel']);
%         digitalInfoChannel
%         auxInfoChannel
        obj.MicrovoltsPerAD = double(obj.electrodeInfoChannel.ConversionFactor); %this is a nx1 array (n channels)
        obj.ZeroADValue=double(obj.electrodeInfoChannel.ADZero);
      
%      Didn't work: obj.channelNumbers = obj.electrodeInfoChannel.ChannelID;
      obj.channelNumbers = 1:length(obj.MicrovoltsPerAD);
      obj.n2s = obj.electrodeInfoChannel.RowIndex;
      obj.channelNames =  obj.electrodeInfoChannel.Label;
%     obj.channelNames = cellfun(@(x) num2str(x), mat2cell(obj.channelNumbers,1,ones(1,numel(obj.channelNumbers))),'UniformOutput',0); 
      obj.analogChannelNumbers=1; %only one non-electrode analog channel - the aux stream
      
      obj.sample_ms = double(obj.electrodeInfoChannel.Tick(1))/1000;
      
      obj.samplingFrequency = 1e3/obj.sample_ms; %Assuming all chanels are the same

%         disp('Extracting time stamp information...');
        obj.timestamps = double(h5read(obj.fullFilename, [obj.pathToRawDataStreamGroup '/ChannelDataTimeStamps']));
%         disp('... done');

%       obj.globalStartTime = double(h5readatt(obj.fullFilename, obj.recNameHD5{1}, 'start_time'));

%       try
%         obj.recordingDuration_ms = double(h5readatt(obj.fullFilename,[obj.recNameHD5{1} '/data'], 'valid_samples'));
%         obj.recordingDuration_ms = 1000 * obj.recordingDuration_ms ./ obj.samplingFrequency;
%         obj.recordingDuration_ms = max(obj.recordingDuration_ms);
%       catch
%         try
%           obj.recordingDuration_ms = double(h5readatt(obj.fullFilename, '/', 'recordingDuration'));
%         catch
%           obj.recordingDuration_ms = 0;
%         end
%       end
      obj.recordingDuration_ms=double(h5readatt(obj.fullFilename,obj.pathToRecording,'Duration'))/1000;
      %%%OR: calculate this from number of sample in channel data, add check by timestamps (number of ticks)
      
      
      obj.unit=obj.electrodeInfoChannel.Unit(1);
      obj.exponent=obj.electrodeInfoChannel.Exponent; 
      rawDataType=obj.electrodeInfoChannel.RawDataType(1);
      databit=obj.electrodeInfoChannel.ADCBits(1);
      if databit==24, databit=32; end %MCs ADC works in 24 bits quantization, but stores in 32 bit format (Try to get aroung this)
      obj.datatype = [char(rawDataType) '32'];                
                
      disp('saving meta data');
      obj.saveMetaData;
    end
  
  end
  
end

function obj = buildN2S(obj)
   obj.n2s=zeros(obj.totalChannels,1);
   for k=1:obj.totalChannels
    channelKLabelIndex=find(chLayoutNumbers==k);   
    obj.n2s(k)=find(strcmp(electrodeInfoChannel.Label,obj.chLayoutNames{channelKLabelIndex}));
   end
end