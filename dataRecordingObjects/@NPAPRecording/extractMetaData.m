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

    %13. ChannelMap

    if isfield(metaAP,'snsGeomMap')
        [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = geomMapToGeom(metaAP);
    elseif isfield(metaAP,'snsShankMap')
        [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = shankMapToGeom(metaAP);
    end

    nchans = numel(xCoord);
    chanPos = zeros(2,nchans);
    for i = 1:nchans
        currX = shankInd(i)*shankPitch + xCoord(i);
        %fprintf( fid, '%d\t%d\t%d\t%d\n', i-1, currX, yCoord(i), shankInd(i));

        chanPos(1,i) = currX;
        chanPos(2,i) = yCoord(i);
    end

    obj.chLayoutPositions = chanPos;


 

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


    % =========================================================
    % Parse snsGeomMap for XY coordinates
    %
    function [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = geomMapToGeom(meta)

        C = textscan(meta.snsGeomMap, '(%d:%d:%d:%d', ...
            'EndOfLine', ')', 'HeaderLines', 1 );
        shankInd = double(cell2mat(C(1)));
        xCoord = double(cell2mat(C(2)));
        yCoord = double(cell2mat(C(3)));
        connected = double(cell2mat(C(4)));

        % parse header for number of shanks
        geomStr = meta.snsGeomMap;
        headStr = extractBefore(geomStr,')(');
        headParts = split(headStr,',');
        nShank = str2double(headParts{2});
        shankWidth = str2double(headParts{4});
        shankPitch = str2double(headParts{3});
    end % geomMapToGeom

% =========================================================
% Get XY coordinates from snsShankMap plus hard coded geom values
%
    function [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = shankMapToGeom(meta)
        % get number of saved AP channels (some early metadata files have a
        % SYNC entry in the snsChanMap
        [nchan,~,~] = ChannelCountsIM(meta);

        C = textscan(meta.snsShankMap, '(%d:%d:%d:%d', ...
            'EndOfLine', ')', 'HeaderLines', 1 );
        shankInd = double(cell2mat(C(1)));
        colInd = double(cell2mat(C(2)));
        rowInd = double(cell2mat(C(3)));
        connected = double(cell2mat(C(4)));

        % trim these to the number of saved channels
        shankInd = shankInd(1:nchan);
        colInd = colInd(1:nchan);
        rowInd = rowInd(1:nchan);
        connected = connected(1:nchan);

        geom = getGeomParams(meta);

        oddRows = logical(mod(rowInd,2));
        evenRows = ~oddRows;
        xCoord = colInd*geom.horzPitch;
        xCoord(evenRows) = xCoord(evenRows) + geom.even_xOff ;
        xCoord(oddRows) = xCoord(oddRows) + geom.odd_xOff;
        yCoord = rowInd*geom.vertPitch;

        nShank = geom.nShank;
        shankWidth = geom.shankWidth;
        shankPitch = geom.shankPitch;
    end % shankMapToGeom

% =========================================================
% Return counts of each imec channel type that compose
% the timepoints stored in binary file.
%
    function [AP,LF,SY] = ChannelCountsIM(meta)
        M = str2num(meta.snsApLfSy);
        AP = M(1);
        LF = M(2);
        SY = M(3);
    end % ChannelCountsIM


% =========================================================
% Return geometry paramters for supported probe types
% These are used to calculate positions from metadata
% that includes only ~snsShankMap
%
    function geom = getGeomParams(meta)
        % create map
        geomTypeMap = makeTypeMap();

        % get probe part number; if absent, this is a 3A
        if isfield(meta,'imDatPrb_pn')
            pn = meta.imDatPrb_pn;
        else
            pn = '3A';
        end

        if geomTypeMap.isKey(pn)
            geomType = geomTypeMap(pn);
        else
            fprintf('unsupported probe part number\n');
            return;
        end

        switch geomType
            case 'np1_stag_70um'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 27;
                geom.odd_xOff = 11;
                geom.horzPitch = 32;
                geom.vertPitch = 20;
                geom.rowsPerShank = 480;
                geom.elecPerShank = 960;
            case 'nhp_lin_70um'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 27;
                geom.odd_xOff = 27;
                geom.horzPitch = 32;
                geom.vertPitch = 20;
                geom.rowsPerShank = 480;
                geom.elecPerShank = 960;
            case 'nhp_stag_125um_med'
                geom.nShank = 1;
                geom.shankWidth = 125;
                geom.shankPitch = 0;
                geom.even_xOff = 27;
                geom.odd_xOff = 11;
                geom.horzPitch = 87;
                geom.vertPitch = 20;
                geom.rowsPerShank = 1368;
                geom.elecPerShank = 2496;
            case 'nhp_stag_125um_long'
                geom.nShank = 1;
                geom.shankWidth = 125;
                geom.shankPitch = 0;
                geom.even_xOff = 27;
                geom.odd_xOff = 11;
                geom.horzPitch = 87;
                geom.vertPitch = 20;
                geom.rowsPerShank = 2208;
                geom.elecPerShank = 4416;
            case 'nhp_lin_125um_med'
                geom.nShank = 1;
                geom.shankWidth = 125;
                geom.shankPitch = 0;
                geom.even_xOff = 11;
                geom.odd_xOff = 11;
                geom.horzPitch = 103;
                geom.vertPitch = 20;
                geom.rowsPerShank = 1368;
                geom.elecPerShank = 2496;
            case 'nhp_lin_125um_long'
                geom.nShank = 1;
                geom.shankWidth = 125;
                geom.shankPitch = 0;
                geom.even_xOff = 11;
                geom.odd_xOff = 11;
                geom.horzPitch = 103;
                geom.vertPitch = 20;
                geom.rowsPerShank = 2208;
                geom.elecPerShank = 4416;
            case 'uhd_8col_1bank'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 14;
                geom.odd_xOff = 14;
                geom.horzPitch = 6;
                geom.vertPitch = 6;
                geom.rowsPerShank = 48;
                geom.elecPerShank = 384;
            case 'uhd_8col_16bank'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 14;
                geom.odd_xOff = 14;
                geom.horzPitch = 6;
                geom.vertPitch = 6;
                geom.rowsPerShank = 768;
                geom.elecPerShank = 6144;
            case 'np2_ss'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 27;
                geom.odd_xOff = 27;
                geom.horzPitch = 32;
                geom.vertPitch = 15;
                geom.rowsPerShank = 640;
                geom.elecPerShank = 1280;
            case 'np2_4s'
                geom.nShank = 4;
                geom.shankWidth = 70;
                geom.shankPitch = 250;
                geom.even_xOff = 27;
                geom.odd_xOff = 27;
                geom.horzPitch = 32;
                geom.vertPitch = 15;
                geom.rowsPerShank = 640;
                geom.elecPerShank = 1280;
            case 'NP1120'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 6.75;
                geom.odd_xOff = 6.75;
                geom.horzPitch = 4.5;
                geom.vertPitch = 4.5;
                geom.rowsPerShank = 192;
                geom.elecPerShank = 384;
            case 'NP1121'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 6.25;
                geom.odd_xOff = 6.25;
                geom.horzPitch = 3;
                geom.vertPitch = 3;
                geom.rowsPerShank = 384;
                geom.elecPerShank = 384;
            case 'NP1122'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 12.5;
                geom.odd_xOff = 12.5;
                geom.horzPitch = 3;
                geom.vertPitch = 3;
                geom.rowsPerShank = 24;
                geom.elecPerShank = 384;
            case 'NP1123'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 10.25;
                geom.odd_xOff = 10.25;
                geom.horzPitch = 4.5;
                geom.vertPitch = 4.5;
                geom.rowsPerShank = 32;
                geom.elecPerShank = 384;
            case 'NP1300'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 11;
                geom.odd_xOff = 11;
                geom.horzPitch = 48;
                geom.vertPitch = 20;
                geom.rowsPerShank = 480;
                geom.elecPerShank = 960;
            case 'NP1200'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 27;
                geom.odd_xOff = 11;
                geom.horzPitch = 32;
                geom.vertPitch = 20;
                geom.rowsPerShank = 64;
                geom.elecPerShank = 128;
            case 'NXT3000'
                geom.nShank = 1;
                geom.shankWidth = 70;
                geom.shankPitch = 0;
                geom.even_xOff = 53;
                geom.odd_xOff = 53;
                geom.horzPitch = 0;
                geom.vertPitch = 15;
                geom.rowsPerShank = 128;
                geom.elecPerShank = 128;
            otherwise
                % shouldn't see this case
                fprintf('unsupported probe part number\n');
                return;
        end
    end %end  getGeomParam

    function M = makeTypeMap()
        % many part numbers have the same geometry parameters ;
        % make a map that pairs geometry type (value) with probe part number (key)
        M = containers.Map('KeyType','char','ValueType','char');

        M('3A') = 'np1_stag_70um';
        M('PRB_1_4_0480_1') = 'np1_stag_70um';
        M('PRB_1_4_0480_1_C') = 'np1_stag_70um';
        M('NP1010') = 'np1_stag_70um';
        M('NP1011') = 'np1_stag_70um';
        M('NP1012') = 'np1_stag_70um';
        M('NP1013') = 'np1_stag_70um';

        M('NP1015') = 'nhp_lin_70um';
        M('NP1015') = 'nhp_lin_70um';
        M('NP1016') = 'nhp_lin_70um';
        M('NP1017') = 'nhp_lin_70um';

        M('NP1020') = 'nhp_stag_125um_med';
        M('NP1021') = 'nhp_stag_125um_med';
        M('NP1030') = 'nhp_stag_125um_long';
        M('NP1031') = 'nhp_stag_125um_long';

        M('NP1022') = 'nhp_lin_125um_med';
        M('NP1032') = 'nhp_lin_125um_long';

        M('NP1100') = 'uhd_8col_1bank';
        M('NP1110') = 'uhd_8col_16bank';

        M('PRB2_1_2_0640_0') = 'np2_ss';
        M('PRB2_1_4_0480_1') = 'np2_ss';
        M('NP2000') = 'np2_ss';
        M('NP2003') = 'np2_ss';
        M('NP2004') = 'np2_ss';

        M('PRB2_4_2_0640_0') = 'np2_4s';
        M('PRB2_4_4_0480_1') = 'np2_4s';
        M('NP2010') = 'np2_4s';
        M('NP2013') = 'np2_4s';
        M('NP2014') = 'np2_4s';

        M('NP1120') = 'NP1120';
        M('NP1121') = 'NP1121';
        M('NP1122') = 'NP1122';
        M('NP1123') = 'NP1123';
        M('NP1300') = 'NP1300';

        M('NP1200') = 'NP1200';
        M('NXT3000') = 'NXT3000';
    end % makeTypeMap

  

end % ReadMeta