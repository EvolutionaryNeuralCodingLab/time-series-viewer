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
    out=regexp(metaAP.fileName,'/','split');
    exp = erase(string(out(end)),'_g0_tcat.imec0.ap.bin');
    obj.recordingName = char(exp); 
    
%     %2. recordingDir % (String) Full directory containing the recorded session
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
    
    obj.chLayoutNumbers = metaAP.snsShankMap; %(MxN) The layout of the channel numbers in physical space arranged in an M by N grid

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
    obj.layoutName = '3B1 staggered';

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
    
    %15. MicrovoltsPerADimec % the digital to analog conversion value
    convImec = Int2Volts(metaAP);
    obj.MicrovoltsPerADimec = convImec;

    %16. MicrovoltsPerADnidq % the digital to analog conversion value
    convNidq = Int2Volts(metaNI);
    obj.MicrovoltsPerADnidq = convNidq;

    %16.0 
    obj.MicrovoltsPerAD = convImec;
    %14.1. 
    obj.convertData2Double = 1; % if data should be converted to double from the original quantization
    
    %-->% Don't know were this info is: ZeroADValue % the digital zero value
    obj.ZeroADValue = 0;


    obj.datatype = 'binary SpikeGLX';        % class of data in the recording

    fprintf('Extracting meta data from: %s...\n',path_r);

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