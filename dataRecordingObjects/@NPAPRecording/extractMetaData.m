%example path: 
% path = '\\132.66.45.127\data\Simon\Anesthesia_exp\Experiment_29_09_22\Insertion1\catgt_Experiment_29_09_22_1_g0';

function obj=extractMetaData(obj)
    
    %Function to get NP metadata
    path_r = obj.recordingDir;

    function [meta] = ReadMetaNP(binName, path_r)

        % Create the matching metafile name
        [dumPath,name,dumExt] = fileparts(binName);
        metaName = strcat(name, '.meta');

        % Parse ini file into cell entries C{1}{i} = C{2}{i}
        fid = fopen(fullfile(path_r, metaName), 'r');
        % -------------------------------------------------------------
        %    Need 'BufSize' adjustment for MATLAB earlier than 2014
        %    C = textscan(fid, '%[^=] = %[^\r\n]', 'BufSize', 32768);
        C = textscan(fid, '%[^=] = %[^\r\n]');
        % -------------------------------------------------------------
        fclose(fid);

        % New empty struct
        meta = struct();

        % Convert each cell entry into a struct entry
        for i = 1:length(C{1})
            tag = C{1}{i};
            if tag(1) == '~'
                % remake tag excluding first character
                tag = sprintf('%s', tag(2:end));
            end
            meta = setfield(meta, tag, C{2}{i});
        end

        

    end % ReadMeta

if iscell(path_r)
    path_r = cell2mat(path_r);
end


    file = dir (path_r);
    filenames = {file.name};
    binFiles = filenames(contains(filenames,".bin"));
    
    %binName, path

    for stream = 1:length(binFiles)
        
        binName = string(binFiles(stream));
        
        meta = ReadMetaNP(binName,path_r);
        
        if contains(binName, '.ap.')
            metaAP = meta;
            metaAP.fileName = binName;
        elseif contains(binName, '.lf.')
            metaLF = meta;
            metaLF.fileName = binName;
        else
            metaNI = meta;
            metaNI.fileName = binName;
        end

    end
    
    % Fill in properties
    % 1. recordingName %(String) The name of the recording
    out=regexp(metaAP.fileName,filesep,'split');
    exp = erase(string(out(end)),'_g0_tcat.imec0.ap.bin');
    obj.recordingName = char(exp); 
    
%   %2. recordingDir % (String) Full directory containing the recorded session
%     out2=regexp(metaAP.fileName,"/"+exp,'split');
%     direc = string(out2(1));
%     obj.recordingDir = direc;
    
    direc = path_r;
    %3. dataFileNames % (Cell 1 x N)  array of N recording data file names
    file = dir (direc);
    filenames = {file.name};
    binFiles = filenames(contains(filenames,".bin"));
    obj.dataFileNames = binFiles;
    
    %4. startDate %(1x1) Start date (time) of Recording (matlab date number format)
    if isfield(metaAP, 'fileCreateTime_original')
        stTime = datetime(metaAP.fileCreateTime_original);
    else
        stTime = datetime(metaAP.fileCreateTime);
    end
    obj.startDate = stTime;

    %5. endDate %(1x1) End date (time) of Recording (matlab long format)
    
    endTime = stTime + seconds(str2double(metaAP.fileTimeSecs));
    obj.endDate = endTime;
    

    %6. samplingFrequencyAP %(1xN) Sampling rate [Hz] of AP band 
    srateAP = str2double(metaAP.imSampRate);
    obj.samplingFrequencyAP = srateAP;

    %7. samplingFrequencyNI %(1xN) Sampling rate [Hz] of NIDQ
    srateNI = str2double(metaNI.niSampRate);
    obj.samplingFrequencyNI = srateNI;
    obj.samplingFrequencyAnalog = srateNI;

    %8. samplingFrequencyLF %(1xN) Sampling rate [Hz] of LF band
    if exist("metaLF","var") 
        srateLF = str2double(metaLF.imSampRate);
        obj.samplingFrequencyLF = srateLF;
    else
        obj.samplingFrequencyLF = 0;
    end

    %6.0 samplingFrequency
    obj.samplingFrequency = srateAP;

    %9. recordingDuration_ms %(1x1) the total duration of the recording in [ms]
    durMS = 1000*str2double(metaAP.fileTimeSecs);
    obj.recordingDuration_ms = durMS;

    
    %10. channelNumbers % (1xN) an array with integer channel numbers (>=1 integer)
    if strcmp(metaAP.snsSaveChanSubset, 'all')
        chans = (1:str2double(metaAP.nSavedChans));
    else
        chans = str2num(metaAP.snsSaveChanSubset);
        chans = chans + 1;
    end
    obj.channelNumbers = 1:numel(chans);

    %10.0. channelNames % (Cell 1xN) a cell array with the N names of the channels 

    obj.channelNames= cellfun(@(x) char(string(x)),mat2cell(chans',ones(1,numel(chans))),'UniformOutput',0);

    %triggerNames %the names of trigger channels (not critical)
    
    %11. analogChannelNumbers % (1xN) the numbers of channels containing non-electrode analog inputs
    if strcmp(metaNI.snsSaveChanSubset, 'all')
        chans1 = (1:str2double(metaNI.nSavedChans));
    else
        chans1 = str2num(metaNI.snsSaveChanSubset);
        chans1 = chans1 + 1;
    end
    obj.analogChannelNumbers = chans1;

    %12. analogChannelNames % (1xN) the numbers of channels containing non-electrode analog inputs
    
    obj.analogChannelNames=cellfun(@(x) ['XA' char(string(x))],mat2cell(chans1(1:end-1)',ones(1,numel(chans1(1:end-1)))),'UniformOutput',0);
    %chans1names(9) = "XD0";

    
    %13.0 
    
    obj.chLayoutNumbers = metaAP.snsChanMap; %(MxN) The layout of the channel numbers in physical space arranged in an M by N grid

    obj.chLayoutNames = (mat2cell(chans,1)); %(Cell MxN)The layout of the channel names in physical space arranged in an M by N grid
    
    % electrodePitch % distance between electrodes (not critical)
    
    %13. chLayoutPositions % (1xN or 2xN or 3xN) array of electrode position in [x or x,y or x,y,z]
    
    ePosy = 20:20:3840;

    ePosY=repelem(ePosy,2);
    
    ePosx =  [43 11 59 27];

    ePosX = repmat(ePosx,1,384/4);

    coor = [ePosX;ePosY];
    obj.chLayoutPositions = coor;

    %14. layoutName %the name of the channel layout (electrode type)
    obj.layoutName = 'staggered-3B1';

    % FUNCTIONS TO RETURN CONVERSION FACTORS:

    %Function to get conversion value:
    % =========================================================
        % Return a multiplicative factor for converting 16-bit
        % file data to voltage. This does not take gain into
        % account. The full conversion with gain is:
        %
        %   dataVolts = dataInt * fI2V / gain.
        %
        % Note that each channel may have its own gain.
        %
    function fI2V = Int2Volts(meta)
        if strcmp(meta.typeThis, 'imec')
            if isfield(meta,'imMaxInt')
                maxInt = str2num(meta.imMaxInt);
            else
                maxInt = 512;
            end
            fI2V = str2double(meta.imAiRangeMax) / maxInt;
        else
            fI2V = str2double(meta.niAiRangeMax) / 32768;
        end
    end % Int2Volts
    
% Return gain arrays for imec channels.
%
% Index into these with original (acquired) channel IDs.
%
    function [APgain,LFgain] = ChanGainsIM(meta)

        if isfield(meta,'imDatPrb_type')
            probeType = str2num(meta.imDatPrb_type);
        else
            probeType = 0;
        end
        if (probeType == 21) || (probeType == 24)
            [AP,LF,~] = ChannelCountsIM(meta);
            % NP 2.0; APgain = 80 for all channels
            APgain = zeros(AP,1,'double');
            APgain = APgain + 80;
            % No LF channels, set gain = 0
            LFgain = zeros(LF,1,'double');
        else
            % 3A or 3B data?
            % 3A metadata has field "typeEnabled" which was replaced
            % with "typeImEnabled" and "typeNiEnabled" in 3B.
            % The 3B imro table has an additional field for the
            % high pass filter enabled/disabled
            if isfield(meta,'typeEnabled')
                % 3A data
                C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d', ...
                    'EndOfLine', ')', 'HeaderLines', 1 );
            else
                % 3B data
                C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d %*s', ...
                    'EndOfLine', ')', 'HeaderLines', 1 );
            end
            APgain = double(cell2mat(C(1)));
            LFgain = double(cell2mat(C(2)));
        end
    end % ChanGainsIM

% =========================================================
% Return array of original channel IDs. As an example,
% suppose we want the imec gain for the ith channel stored
% in the binary data. A gain array can be obtained using
% ChanGainsIM() but we need an original channel index to
% do the look-up. Because you can selectively save channels
% the ith channel in the file isn't necessarily the ith
% acquired channel, so use this function to convert from
% ith stored to original index.
%
% Note: In SpikeGLX channels are 0-based, but MATLAB uses
% 1-based indexing, so we add 1 to the original IDs here.
%
    function chans = OriginalChans(meta)
        if strcmp(meta.snsSaveChanSubset, 'all')
            chans = (1:str2double(meta.nSavedChans));
        else
            chans = str2num(meta.snsSaveChanSubset);
            chans = chans + 1;
        end
    end % OriginalChans

    [APgain,LFgain] = ChanGainsIM(metaAP);
    APgain= APgain(1); %takes gain of first channel. Asumes all channels have the same gain).

    %15. MicrovoltsPerADimec % the digital to analog conversion value AP
    obj.MicrovoltsPerAD = repmat(Int2Volts(metaAP)/APgain,1,length(chans));

    %15.2 MicrovoltsPerADimec % the digital to analog conversion value LF
    LFgain= LFgain(1);
    obj.MicrovoltsPerAD_LF = repmat(Int2Volts(metaAP)/LFgain,1,length(chans));

    %Functions required to get gain from dinaq channels
        % =========================================================
    % Return gain for ith channel stored in the nidq file.
    %
    % ichan is a saved channel index, rather than an original
    % (acquired) index.
    %
    function gain = ChanGainNI(ichan, savedMN, savedMA, meta)
        if ichan <= savedMN
            gain = str2double(meta.niMNGain);
        elseif ichan <= savedMN + savedMA
            gain = str2double(meta.niMAGain);
        else
            gain = 1;
        end
    end % ChanGainNI

    % =========================================================
    % Return counts of each nidq channel type that compose
    % the timepoints stored in binary file.
    %
    function [MN,MA,XA,DW] = ChannelCountsNI(meta)
        M = str2num(meta.snsMnMaXaDw);
        MN = M(1);
        MA = M(2);
        XA = M(3);
        DW = M(4);
    end % ChannelCountsNI

    %16. MicrovoltsPerADnidq % the digital to analog conversion value
     [MN,MA] = ChannelCountsNI(metaNI);
     fI2V = Int2Volts(metaNI);

      for i = 1:length(chans1)   % index into timepoint
            NIconv(i) = fI2V / ChanGainNI(i, MN, MA, meta);
  
      end

    convNidq = Int2Volts(metaNI)/NIconv(1); %takes gain of first channel. Asumes all channels have the same gain).
    obj.MicrovoltsPerADAnalog= convNidq; 

    
    
    %17. Zero value for binary recordings
    obj.ZeroADValue = 0;
    obj.ZeroADValueAnalog = 0;

    %18. Object type
    obj.datatype = 'int16';        % class of data in the recording

    fprintf('Extracting meta data from: %s...\n',path_r);


    %Channels 
    obj.nSavedChansNidq = metaNI.nSavedChans;
    obj.nSavedChansImec = metaAP.nSavedChans;

    obj.fileSizeBytesNI = metaNI.fileSizeBytes;
    obj.fileSizeBytesAP = metaAP.fileSizeBytes;
    
    if exist("metaLF","var") 
        obj.fileSizeBytesLF = metaLF.fileSizeBytes;
    else
        obj.fileSizeBytesLF = 0;
    end
            
    %disp('saving meta data');
    %obj.saveMetaData(metaAP);

end % ReadMeta