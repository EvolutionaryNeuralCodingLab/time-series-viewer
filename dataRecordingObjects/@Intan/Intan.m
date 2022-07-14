classdef Intan < dataRecording
    properties
        %{
        Properties of parent class: dataRecording
        recordingDir %Full directory containing the recorded session
        startDate %Start date (time) of Recording (matlab long format)
        endDate %End date (time) of Recording (matlab long format)
        samplingFrequency %Sampling rate
        dataFileNames %Name of recording data files
        channelNames %a cell array with the names of the channels
        channelNumbers %an array with integer channel numbers
        triggerNames %the names of trigger channels
        recordingDuration_ms %the total duration of the recording
        %}
    end
    properties (Constant)
        defaultLocalDir='C:\Users\Tulip\Documents\Academic\Post-Doc\Experiments'; %Default directory from which search starts
        fileExtension='rhd';
        maxNumberOfDigitalChannels=4;
        
        signalBits = 16; %the quantization of the sampling card
        numberOfCharFromEndToBaseName=7;
    end
    properties (Hidden)
        startDateLocal
        endDateLocal
        recordingDurationLocal_ms
        cumStart
        cumEnd
        samplesPerMS
        data_file_main_version_number
        data_file_secondary_version_number
        num_data_blocks
        bytes_remaining
        bytes_per_block
        num_amplifier_channels
        num_aux_input_channels
        num_supply_voltage_channels
        num_board_adc_channels
        num_board_dig_in_channels
        num_board_dig_out_channels
        num_temp_sensor_channels
        dataBytesStart
        num_amplifier_samples
        num_aux_input_samples
        num_supply_voltage_samples
        num_board_adc_samples;
        num_board_dig_in_samples;
        num_board_dig_out_samples;
        chanElectrodes
        electrNum
        chan2elect
        % layoutName
    end
    methods
        
        
        function [V_uV,t_ms]=getData(obj,channels,startTime_ms,window_ms,name)
            %Extract MCRack recording data from file to memory
            %Usage: [V_uV,T_ms]=obj.getData(channels,startTime_ms,window_ms,name);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the stream (if not entered, default name is used)
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        T_ms - A time vector relative to recording start (t=0 at start)
            
            %transform channel request to actual intan channel number (different for
            %multiple electrodes)
            if size(obj,2)>1
                error('You have to select a single object from the object array!');
            end
            
            %channels=obj.chan2elect(channels);
            
            conversionFactor=1/1000*obj.samplingFrequency;
            endTime_ms=startTime_ms+window_ms; %no need to conversion factor
            windowSamples=floor(window_ms*conversionFactor);
            
            %Intan interleaves data blocks with other things, to load just
            %data we need to skip over the other stuff
            dataBytesPerBlock = 60*2*obj.num_amplifier_channels;
            dataPointsPerBlock =  60;
            nonDataBytesPerBlock =  obj.bytes_per_block-dataBytesPerBlock;
            
            
            nTrials=length(startTime_ms);
            %make it a number of blocks big, load the blocks for all trials and files, trim the beginning and end
            V_uV=zeros(numel(channels),nTrials,windowSamples);
            amplifier_index = 1;
            
            cumStart=[-Inf obj.cumStart obj.cumEnd(end)];
            cumEnd=[0 obj.cumEnd Inf];
            
            % if obj.multifileMode %this mode currently does not support extraction from edges of the recording
            for i=1:nTrials
                tmpStartTime=startTime_ms(i);
                dataTally=1;
                pFileStart=find(startTime_ms(i)>=cumStart,1,'last');
                pFileEnd=find((startTime_ms(i)+window_ms)<=cumEnd,1,'first');
                
                for f=pFileStart:pFileEnd
                    if f>1 && f<=(obj.nRecordings+1) % data in inside recording range
                        
                        tmpEndTime=min([cumEnd(f) endTime_ms(i)]);    %start and end time in the file f
                        
                        if f == pFileStart
                            startSampleLocal=round((tmpStartTime-cumStart(f))*obj.samplesPerMS); %for sample number in neuro-share - like tmpStartTime
                        else
                            startSampleLocal=0;
                        end
                        
                        windowSamplesLocal=((tmpEndTime-cumEnd(f-1))*obj.samplesPerMS)-startSampleLocal;
                        endSampleLocal=startSampleLocal+windowSamplesLocal;
                        
                        filename=[obj.recordingDir filesep obj.dataFileNames{f-1}];
                        fid = fopen(filename,'r');
                        s = dir(filename);
                        filesize = s.bytes;
                        
                        startBlock=floor(startSampleLocal/dataPointsPerBlock);
                        endBlock=ceil(endSampleLocal/dataPointsPerBlock);
                        
                        tmpStartSample=(startBlock)*dataPointsPerBlock;
                        
                        %load the data from startBlock to endBlock. This data
                        %lines up with certain sample times. Trim this to the
                        %desired start and end samples
                        
                        fseek(fid,startBlock*obj.bytes_per_block+obj.dataBytesStart+240,'bof');
                        dataBytes=60*obj.num_amplifier_channels;
                        data = fread(fid,[60*obj.num_amplifier_channels,endBlock-startBlock+1],[int2str(dataBytes) '*uint16'],obj.bytes_per_block-obj.num_amplifier_channels*60*2);
                        
                        data=reshape(data, 60 ,obj.num_amplifier_channels,size(data,2));
                        data=permute(data,[1 3 2]);
                        data=reshape(data,size(data,1)*size(data,2),size(data,3));
                        
                        
                        data=obj.MicrovoltsPerAD * (data - 32768); %scale the data
                        
                        indicesToKeep=(startSampleLocal-tmpStartSample+1):(startSampleLocal-tmpStartSample+windowSamplesLocal);
                        
                        V_uV(:, i, dataTally:(dataTally+numel(indicesToKeep)-1))= data(indicesToKeep,channels)';
                        dataTally=dataTally+size(data(indicesToKeep,:),1);
                        
                        fclose(fid);
                    else
                        V_uV(:, i, dataTally:end) = zeros(size(V_uV,1),1,size(V_uV,3)-dataTally+1);
                    end
                end
            end
            
            t_ms=(1:size(V_uV,3))/obj.samplesPerMS;
            
        end
        
        function [D,T_ms]=getDigitalData(obj,startTime_ms,window_ms,name)
            %Extract Intan digital data from file to memory
            %Usage: [D,T_ms]=getDigitalData(startTime_ms,window_ms,name)
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the stream (if not entered, default name is used)
            %Output: D - A 3D matrix [nChannels x nTrials x nSamples] with digitalData waveforms across specified channels and trials
            %        T_ms - A time vector relative to recording start (t=0 at start)
            
              if size(obj,2)>1
                obj=obj(1);
            end
            
            if nargin==3
                pDigitalEntity=1:obj.num_board_dig_in_channels; %take all triggers in recording
            elseif nargin==4
                error('Entering digial data name is  currently not supported for intan data');
                %this option should be revised because currently all parameters are derived from the raw data stream
            elseif nargin==2
                error('Not enough input variables: window is required when entering time frame');
            elseif nargin==1
                startTime_ms=0;
                window_ms=obj.recordingDuration_ms;
                pDigitalEntity=1:obj.num_board_dig_in_channels; %take all triggers in recording
            else
                error('method getDigitalData was not used correctly: wrong number of inputs');
            end
            
                    
            
            conversionFactor=1/1000*obj.samplingFrequency;
            endTime_ms=startTime_ms+window_ms; %no need to conversion factor
            windowSamples=round(window_ms*conversionFactor);
            
            %Intan interleaves data blocks with other things, to load just
            %data we need to skip over the other stuff
            dataBytesPerBlock = 60*2;
            dataPointsPerBlock =  60;
            nonDataBytesPerBlock =  obj.bytes_per_block-dataBytesPerBlock;
            
            
            nTrials=length(startTime_ms);
            %make it a number of blocks big, load the blocks for all trials and files, trim the beginning and end
            
            D=false(obj.num_board_dig_in_channels,nTrials,windowSamples);
            amplifier_index = 1;
            
            
            cumStart=[-Inf obj.cumStart obj.cumEnd(end)];
            cumEnd=[0 obj.cumEnd Inf];
            
            % if obj.multifileMode %this mode currently does not support extraction from edges of the recording
            for i=1:nTrials
                
                tmpStartTime=startTime_ms(i);
                dataTally=1;
                pFileStart=find(startTime_ms(i)>=cumStart,1,'last');
                pFileEnd=find((startTime_ms(i)+window_ms)<=cumEnd,1,'first');
                
                for f=pFileStart:pFileEnd
                    if f>1 && f<=(obj.nRecordings+1) % data in inside recording range
                        
                        tmpEndTime=min([cumEnd(f) endTime_ms(i)]);    %start and end time in the file f
                        
                        if f == pFileStart
                            startSampleLocal=round((tmpStartTime-cumStart(f))*obj.samplesPerMS); %for sample number in neuro-share - like tmpStartTime
                        else
                            startSampleLocal=0;
                        end
                        
                        windowSamplesLocal=(tmpEndTime-cumEnd(f-1))*obj.samplesPerMS-startSampleLocal;
                        endSampleLocal=startSampleLocal+windowSamplesLocal;
                        
                        filename=[obj.recordingDir filesep obj.dataFileNames{f-1}];
                        fid = fopen(filename,'r');
                        s = dir(filename);
                        filesize = s.bytes;
                        
                        startBlock=floor(startSampleLocal/dataPointsPerBlock);
                        endBlock=ceil(endSampleLocal/dataPointsPerBlock);
                        
                        tmpStartSample=(startBlock)*dataPointsPerBlock;
                        
                        %load the data from startBlock to endBlock. This data
                        %lines up with certain sample times. Trim this to the
                        %desired start and end samples
                        bytesBeforeRelevantShit=60 * 2 * 2+...
                            60 * 2 * obj.num_amplifier_channels+...
                            15 * 2 * obj.num_aux_input_channels+...
                            1 * 2 * obj.num_supply_voltage_channels+...
                            1 * 2 * obj.num_temp_sensor_channels+...
                            60 * 2 * obj.num_board_adc_channels;
                        
                        %previous for loop was replaced by Mark to increase speed
                        fseek(fid,startBlock*obj.bytes_per_block + obj.dataBytesStart + bytesBeforeRelevantShit,'bof');
                        data = fread(fid,[60,endBlock-startBlock+1],'60*uint16',obj.bytes_per_block-2*60);
                        
                        data=reshape(data,1,size(data,1)*size(data,2));
                        
                        %seperate into different channels
                        board_dig_in_data=false(obj.num_board_dig_in_channels,length(data));
                        for ii=1:obj.num_board_dig_in_channels
                            mask = 2^(ii-1) * ones(1,length(data));
                            board_dig_in_data(ii, :) = (bitand(data, mask) > 0);
                        end
                        
                        indicesToKeep=((startSampleLocal+1)-tmpStartSample):((startSampleLocal+1)-tmpStartSample+windowSamplesLocal-1);
                        
                        D(:,i,dataTally:(dataTally+numel(indicesToKeep)-1))= board_dig_in_data(:,indicesToKeep);
                        dataTally=dataTally+size(data(:,indicesToKeep),2);
                        
                    else
                        D( :,i, dataTally:end) = false(size(D,1),size(D,3)-dataTally+1);
                    end
                end
                fclose(fid);
            end
            
            T_ms=(1:size(D,3))/obj.samplesPerMS;
            
        end
        
        function [T_ms]=getTrigger(obj,startTime_ms,window_ms,name)
            %Extract triggers from recording
            %Usage : [T_ms]=obj.getTrigger(,startTime_ms,window_ms,name)
            %Input : startTime_ms - start time [ms].
            %        window_ms - the window duration [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        name
            %Output: T_ms - trigger times [ms] - different triggers are arranged in a cell array
            
                if size(obj,2)>1
                obj=obj(1);
            end
            
            if obj.num_board_dig_in_channels<1
                disp('Warning!!!! No trigger or digital data entities in recording');
                T_ms={};
                return;
            end
                  
            if nargin==4
                disp(['Triggers in bits: ' num2str(name) ' will be extracted']);
                digitalBits=name;
            elseif nargin==3
                digitalBits=1:16;
            elseif nargin==1
                startTime_ms=0;
                window_ms=obj.recordingDuration_ms;
                digitalBits=1:16;
            end
            
            endTime_ms=startTime_ms+window_ms; %no need to conversion factor
            
            pFileStart=find(startTime_ms>=obj.cumStart,1,'last');
            pFileEnd=max([find((startTime_ms+window_ms)<=obj.cumEnd,1,'first'),1]);
            nFiles=pFileEnd-pFileStart+1;
            
            nTriggersDigital=2*numel(digitalBits);
            T=cell(nTriggersDigital,nFiles);
            
            tmpStartTime=startTime_ms;
            for i=pFileStart:pFileEnd
                tmpEndTime=min([obj.cumEnd(i) endTime_ms]);
                
                [D,Ttmp]=getDigitalData(obj,tmpStartTime,tmpEndTime-tmpStartTime);
                validChannels=any(any(D,3),2);
                D=D(validChannels,:,:);
                for j=1:size(D,1)
                    T{2*j-1,i}=Ttmp(find(diff(squeeze(D(j,:,:)))>0))+obj.cumStart(i);
                    T{2*j,i}=Ttmp(find(diff(squeeze(D(j,:,:)))<0))+obj.cumStart(i);
                end
                tmpStartTime=tmpEndTime;
            end
            
            for i=1:nTriggersDigital
                T_ms{i}=cell2mat(T(i,:));
            end
            T_ms(cellfun(@isempty,T_ms))=[];
            
        end
        
    end
    
    %% sam
    methods (Hidden)
        %class constructor
        function obj = Intan(recordingFile,A)
            
            obj.datatype='uint16';
            obj.MicrovoltsPerAD=0.195;
            
            %object array for multiple electrodes
            if nargin == 2
                if A~=1 %to stop infinite recursion!
                    for x=1:A
                        obj(x)=Intan(recordingFile,1);
                        obj(x).signalBits=x;
                        obj(x).channelNumbers=1:numel(find(obj(x).chanElectrodes == obj(x).electrNum(x)));
                        obj(x).channelNames=obj(x).channelNumbers;
                        obj(x).chan2elect=find(obj(x).chanElectrodes == obj(x).electrNum(x));
                        if iscell(obj(x).recordingName)
                        obj(x).recordingName=['elect' int2str(x) '_'  obj(x).recordingName{1}];
                        end
                        if ~isempty(obj(x).layoutName)
                            obj(x).layoutName=obj(x).layoutName{x};
                            obj(x).chLayoutNumbers= obj(x).chLayoutNumbers{x};
                            obj(x).chLayoutNames= obj(x).chLayoutNames{x};
                        end
                        
                    end
                    return
                end
            end
            
            %get data files
            if nargin==0
                recordingFile=[];
            elseif nargin>2
                disp('Object was not constructed since too many parameters were given at construction');
                return;
            end
            obj=obj.getRecordingFiles(recordingFile,obj.fileExtension);
            
            %open first file for reading and getting file data that applies
            %to all files
            for i=1:obj.nRecordings
                
                filename=[obj.recordingDir filesep obj.dataFileNames{i}];
                fid = fopen(filename,'r');
                s = dir(filename);
                filesize = s.bytes;
                
                % Check 'magic number' at beginning of file to make sure this is an Intan
                % Technologies RHD2000 data file.
                magic_number = fread(fid, 1, 'uint32');
                if magic_number ~= hex2dec('c6912702')
                    error('Unrecognized file type.');
                end
                
                % Read version number.
                obj.data_file_main_version_number = fread(fid, 1, 'int16');
                obj.data_file_secondary_version_number = fread(fid, 1, 'int16');
                
                if i==1
                    if nargin==1
                        fprintf(1, 'Reading Intan Technologies RHD2000 Data File, Version %d.%d\n', ...
                            obj.data_file_main_version_number, obj.data_file_secondary_version_number);
                    end
                end
                
                % Read information of sampling rate and amplifier frequency settings.
                sample_rate = fread(fid, 1, 'single');
                dsp_enabled = fread(fid, 1, 'int16');
                actual_dsp_cutoff_frequency = fread(fid, 1, 'single');
                actual_lower_bandwidth = fread(fid, 1, 'single');
                actual_upper_bandwidth = fread(fid, 1, 'single');
                
                desired_dsp_cutoff_frequency = fread(fid, 1, 'single');
                desired_lower_bandwidth = fread(fid, 1, 'single');
                desired_upper_bandwidth = fread(fid, 1, 'single');
                
                
                % This tells us if a software 50/60 Hz notch filter was enabled during
                % the data acquisition.
                notch_filter_mode = fread(fid, 1, 'int16');
                notch_filter_frequency = 0;
                if (notch_filter_mode == 1)
                    notch_filter_frequency = 50;
                elseif (notch_filter_mode == 2)
                    notch_filter_frequency = 60;ss
                end
                
                desired_impedance_test_frequency = fread(fid, 1, 'single');
                actual_impedance_test_frequency = fread(fid, 1, 'single');
                
                
                % Place notes in data strucure
                notes = struct( ...
                    'note1', fread_QString(fid), ...
                    'note2', fread_QString(fid), ...
                    'note3', fread_QString(fid) );
                
                % If data file is from GUI v1.1 or later, see if temperature sensor data
                % was saved.
                obj.num_temp_sensor_channels = 0;
                if ((obj.data_file_main_version_number == 1 && obj.data_file_secondary_version_number >= 1) ...
                        || (obj.data_file_main_version_number > 1))
                    obj.num_temp_sensor_channels = fread(fid, 1, 'int16');
                end
                
                % If data file is from GUI v1.3 or later, load eval board mode.
                eval_board_mode = 0;
                if ((obj.data_file_main_version_number == 1 && obj.data_file_secondary_version_number >= 3) ...
                        || (obj.data_file_main_version_number > 1))
                    eval_board_mode = fread(fid, 1, 'int16');
                end
                
                % Place frequency-related information in data structure.
                frequency_parameters = struct( ...
                    'amplifier_sample_rate', sample_rate, ...
                    'aux_input_sample_rate', sample_rate / 4, ...
                    'supply_voltage_sample_rate', sample_rate / 60, ...
                    'board_adc_sample_rate', sample_rate, ...
                    'board_dig_in_sample_rate', sample_rate, ...
                    'desired_dsp_cutoff_frequency', desired_dsp_cutoff_frequency, ...
                    'actual_dsp_cutoff_frequency', actual_dsp_cutoff_frequency, ...
                    'dsp_enabled', dsp_enabled, ...
                    'desired_lower_bandwidth', desired_lower_bandwidth, ...
                    'actual_lower_bandwidth', actual_lower_bandwidth, ...
                    'desired_upper_bandwidth', desired_upper_bandwidth, ...
                    'actual_upper_bandwidth', actual_upper_bandwidth, ...
                    'notch_filter_frequency', notch_filter_frequency, ...
                    'desired_impedance_test_frequency', desired_impedance_test_frequency, ...
                    'actual_impedance_test_frequency', actual_impedance_test_frequency );
                
                % Define data structure for spike trigger settings.
                spike_trigger_struct = struct( ...
                    'voltage_trigger_mode', {}, ...
                    'voltage_threshold', {}, ...
                    'digital_trigger_channel', {}, ...
                    'digital_edge_polarity', {} );
                
                new_trigger_channel = struct(spike_trigger_struct);
                spike_triggers = struct(spike_trigger_struct);
                
                % Define data structure for data channels.
                channel_struct = struct( ...
                    'native_channel_name', {}, ...
                    'custom_channel_name', {}, ...
                    'native_order', {}, ...
                    'custom_order', {}, ...
                    'board_stream', {}, ...
                    'chip_channel', {}, ...
                    'port_name', {}, ...
                    'port_prefix', {}, ...
                    'port_number', {}, ...
                    'electrode_impedance_magnitude', {}, ...
                    'electrode_impedance_phase', {} );
                
                new_channel = struct(channel_struct);
                
                % Create structure arrays for each type of data channel.
                amplifier_channels = struct(channel_struct);
                aux_input_channels = struct(channel_struct);
                supply_voltage_channels = struct(channel_struct);
                board_adc_channels = struct(channel_struct);
                board_dig_in_channels = struct(channel_struct);
                board_dig_out_channels = struct(channel_struct);
                
                amplifier_index = 1;
                aux_input_index = 1;
                supply_voltage_index = 1;
                board_adc_index = 1;
                board_dig_in_index = 1;
                board_dig_out_index = 1;
                
                % Read signal summary from data file header.
                
                number_of_signal_groups = fread(fid, 1, 'int16');
                
                for signal_group = 1:number_of_signal_groups
                    signal_group_name = fread_QString(fid);
                    signal_group_prefix = fread_QString(fid);
                    signal_group_enabled = fread(fid, 1, 'int16');
                    signal_group_num_channels = fread(fid, 1, 'int16');
                    signal_group_num_amp_channels = fread(fid, 1, 'int16');
                    
                    if (signal_group_num_channels > 0 && signal_group_enabled > 0)
                        new_channel(1).port_name = signal_group_name;
                        new_channel(1).port_prefix = signal_group_prefix;
                        new_channel(1).port_number = signal_group;
                        for signal_channel = 1:signal_group_num_channels
                            new_channel(1).native_channel_name = fread_QString(fid);
                            new_channel(1).custom_channel_name = fread_QString(fid);
                            new_channel(1).native_order = fread(fid, 1, 'int16');
                            new_channel(1).custom_order = fread(fid, 1, 'int16');
                            signal_type = fread(fid, 1, 'int16');
                            channel_enabled = fread(fid, 1, 'int16');
                            new_channel(1).chip_channel = fread(fid, 1, 'int16');
                            new_channel(1).board_stream = fread(fid, 1, 'int16');
                            new_trigger_channel(1).voltage_trigger_mode = fread(fid, 1, 'int16');
                            new_trigger_channel(1).voltage_threshold = fread(fid, 1, 'int16');
                            new_trigger_channel(1).digital_trigger_channel = fread(fid, 1, 'int16');
                            new_trigger_channel(1).digital_edge_polarity = fread(fid, 1, 'int16');
                            new_channel(1).electrode_impedance_magnitude = fread(fid, 1, 'single');
                            new_channel(1).electrode_impedance_phase = fread(fid, 1, 'single');
                            
                            if (channel_enabled)
                                switch (signal_type)
                                    case 0
                                        amplifier_channels(amplifier_index) = new_channel;
                                        spike_triggers(amplifier_index) = new_trigger_channel;
                                        amplifier_index = amplifier_index + 1;
                                    case 1
                                        aux_input_channels(aux_input_index) = new_channel;
                                        aux_input_index = aux_input_index + 1;
                                    case 2
                                        supply_voltage_channels(supply_voltage_index) = new_channel;
                                        supply_voltage_index = supply_voltage_index + 1;
                                    case 3
                                        board_adc_channels(board_adc_index) = new_channel;
                                        board_adc_index = board_adc_index + 1;
                                    case 4
                                        board_dig_in_channels(board_dig_in_index) = new_channel;
                                        board_dig_in_index = board_dig_in_index + 1;
                                    case 5
                                        board_dig_out_channels(board_dig_out_index) = new_channel;
                                        board_dig_out_index = board_dig_out_index + 1;
                                    otherwise
                                        error('Unknown channel type');
                                end
                            end
                            
                        end
                    end
                end
                
                obj.num_amplifier_channels = amplifier_index - 1;
                obj.num_aux_input_channels = aux_input_index - 1;
                obj.num_supply_voltage_channels = supply_voltage_index - 1;
                obj.num_board_adc_channels = board_adc_index - 1;
                obj.num_board_dig_in_channels = board_dig_in_index - 1;
                obj.num_board_dig_out_channels = board_dig_out_index - 1;
                
                
                if i==1
                    if nargin==1
                        fprintf(1, 'Found %d amplifier channel%s.\n', ...
                            obj.num_amplifier_channels, plural( obj.num_amplifier_channels));
                        fprintf(1, 'Found %d auxiliary input channel%s.\n', ...
                            obj.num_aux_input_channels, plural( obj.num_aux_input_channels));
                        fprintf(1, 'Found %d supply voltage channel%s.\n', ...
                            obj.num_supply_voltage_channels, plural( obj.num_supply_voltage_channels));
                        fprintf(1, 'Found %d board ADC channel%s.\n', ...
                            obj.num_board_adc_channels, plural( obj.num_board_adc_channels));
                        fprintf(1, 'Found %d board digital input channel%s.\n', ...
                            obj.num_board_dig_in_channels, plural( obj.num_board_dig_in_channels));
                        fprintf(1, 'Found %d board digital output channel%s.\n', ...
                            obj.num_board_dig_out_channels, plural( obj.num_board_dig_out_channels));
                        fprintf(1, 'Found %d temperature sensors channel%s.\n', ...
                            obj.num_temp_sensor_channels, plural( obj.num_temp_sensor_channels));
                        fprintf(1, '\n');
                    end
                end
                
                % Each data block contains 60 amplifier samples.
                obj.bytes_per_block = 60 * 4;  % timestamp data
                obj.bytes_per_block = obj.bytes_per_block + 60 * 2 * obj.num_amplifier_channels;
                % Auxiliary inputs are sampled 4x slower than amplifiers
                obj.bytes_per_block = obj.bytes_per_block + 15 * 2 * obj.num_aux_input_channels;
                % Supply voltage is sampled 60x slower than amplifiers
                obj.bytes_per_block = obj.bytes_per_block + 1 * 2 * obj.num_supply_voltage_channels;
                % Board analog inputs are sampled at same rate as amplifiers
                obj.bytes_per_block = obj.bytes_per_block + 60 * 2 * obj.num_board_adc_channels;
                % Board digital inputs are sampled at same rate as amplifiers
                if (obj.num_board_dig_in_channels > 0)
                    obj.bytes_per_block = obj.bytes_per_block + 60 * 2;
                end
                % Board digital outputs are sampled at same rate as amplifiers
                if (obj.num_board_dig_out_channels > 0)
                    obj.bytes_per_block = obj.bytes_per_block + 60 * 2;
                end
                % Temp sensor is sampled 60x slower than amplifiers
                if (obj.num_temp_sensor_channels > 0)
                    obj.bytes_per_block = obj.bytes_per_block + 1 * 2 * obj.num_temp_sensor_channels;
                end
                
                % How many data blocks remain in this file?
                obj.bytes_remaining(i) = filesize - ftell(fid);
                
                obj.num_data_blocks(i) = obj.bytes_remaining(i) / obj.bytes_per_block;
                
                obj.num_amplifier_samples(i) = 60 * obj.num_data_blocks(i);
                obj.num_aux_input_samples(i) = 15 * obj.num_data_blocks(i);
                obj.num_supply_voltage_samples(i) = 1 * obj.num_data_blocks(i);
                obj.num_board_adc_samples(i) = 60 * obj.num_data_blocks(i);
                obj.num_board_dig_in_samples(i) = 60 * obj.num_data_blocks(i);
                obj.num_board_dig_out_samples(i) = 60 * obj.num_data_blocks(i);
                
                
                record_time = obj.num_amplifier_samples(i) / sample_rate;
                
                
                obj.recordingDurationLocal_ms(i)=record_time*1000;
                
                %get date
                if int2str(sum(obj.dataFileNames{i}=='_'))~='0'
                    breaks=find(obj.dataFileNames{i}=='_');
                    ymd=obj.dataFileNames{i}((breaks(1)+1):(breaks(2)-1));
                    hms=obj.dataFileNames{i}((breaks(2)+1):(breaks(2)+6));
                end
                obj.startDateLocal(i)=datenum([str2num(ymd(1:2)),str2num(ymd(3:4)),str2num(ymd(5:6)), str2num(hms(1:2)), str2num(hms(3:4)), str2num(hms(5:6))]);
                obj.endDateLocal(i)=obj.startDateLocal(i)+datenum(0,0,0,0,0,obj.recordingDurationLocal_ms(i)/1000);
                
            end
            
            obj.recordingDuration_ms=sum(obj.recordingDurationLocal_ms);
            obj.cumEnd=cumsum(obj.recordingDurationLocal_ms);
            obj.cumStart=[0 obj.cumEnd(1:end-1)];
            obj.startDate= obj.startDateLocal(1);
            obj.endDate=obj.endDateLocal(end);
            
            if obj.nRecordings>1
                obj.multifileMode=true;
            else
                obj.multifileMode=false;
            end
            %check that files were entered in the correct order
            if obj.multifileMode
                
                if min(diff(obj.startDateLocal))<0 %changed it a bit, not sure what Mark's problem was with this line
                    disp('Error: files were not entered in their chronological order!!! Please reorder before running');
                    return;
                end
                
            end
            
            % These properties are considered to be the same for all recorded files and extracted only from the first file
            obj.samplingFrequency=sample_rate; %[Hz]
            obj.samplesPerMS=obj.samplingFrequency/1000;
            
            % which channel comes from which electros
            for chan=1:numel(amplifier_channels)
                obj.chanElectrodes(chan)=amplifier_channels(chan).port_number;
            end
            
            %load layout
            obj=obj.loadChLayout;
            
            if numel(~isnan(obj.chLayoutNumbers))> obj.chanElectrodes(1)
                obj.channelNames={amplifier_channels.custom_channel_name};
                obj.channelNumbers=1:numel(amplifier_channels);
                singleElectrode=true;
            else
                obj.channelNumbers=cell2mat({amplifier_channels.custom_order})+1;
                obj.channelNames=mat2cell(obj.channelNumbers,1,ones(1,numel(obj.channelNumbers)));
            end
            obj.dataBytesStart = ftell(fid);
            
            
            % functionality for multiple electrodes
            for chan=1:numel(amplifier_channels)
                obj.chanElectrodes(chan)=amplifier_channels(chan).port_number;
            end
            
            obj.electrNum=unique(obj.chanElectrodes);
            obj.chan2elect=find(obj.chanElectrodes == obj.electrNum(1));
            
            
            fclose(fid);
            
            if nargin==1
                if numel(obj.electrNum)>1 & ~singleElectrode
                    A=numel(obj.electrNum);
                    obj=Intan(recordingFile,A);
                end
            end
            
        end
    end
end