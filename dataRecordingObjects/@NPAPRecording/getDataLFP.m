function [V_uV,t_ms]=getDataLFP(obj,channels,startTime_ms,window_ms) 

%obj,channels,startTime_ms,window_ms


            %Extract recording data from file to memory
            %Usage: [V_uV,t_ms]=obj.getData(channels,startTime_ms,window_ms);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the recording (if empty takes the default name)
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        t_ms - A time vector relative to recording start (t=0 at start)


            %%%%% IMPLEMENTATION:
    binName = string(obj.dataFileNames(contains(obj.dataFileNames,'lf')));
    
    path = obj.recordingDir;
    
    meta = ReadMeta(binName, path);
    
    
    chanList = channels;
    
    nSamp = round((window_ms/1000)*str2double(meta.imSampRate));
    
    V_uV = zeros(length(channels),length(startTime_ms),nSamp);
    
    for trials = 1:length(startTime_ms)

        t0= startTime_ms(trials)/1000;
        samp0 = round(t0*str2double(meta.imSampRate));
    
        dataArray = ReadBin(samp0, nSamp, meta, binName, path);
        
        dataArray = dataArray(chanList,:);

        V_uV(:,trials,:) = GainCorrectIM(dataArray*1000000, chanList, meta); %convert to microvolts
    end
    
    
%     nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
%     
%     samples =samp0:nSamp+samp0;
%     
%     t_ms = vpa((samples/str2double(meta.imSampRate))*1000);

   t_ms = (1:nSamp)/str2double(meta.imSampRate)*1000;



% =========================
% General Utility Functions
% =========================


% =========================================================
% Parse ini file returning a structure whose field names
% are the metadata left-hand-side tags, and whose right-
% hand-side values are MATLAB strings. We remove any
% leading '~' characters from tags because MATLAB uses
% '~' as an operator.
%
% If you're unfamiliar with structures, the benefit
% is that after calling this function you can refer
% to metafile items by name. For example:
%
%   meta.fileCreateTime  // file create date and time
%   meta.nSavedChans     // channels per timepoint
%
% All of the values are MATLAB strings, but you can
% obtain a numeric value using str2double(meta.nSavedChans).
% More complicated parsing of values is demonstrated in the
% utility functions below.
%
    function [meta] = ReadMeta(binName, path)

        % Create the matching metafile name
        [dumPath,name,dumExt] = fileparts(binName);
        metaName = strcat(name, '.meta');

        % Parse ini file into cell entries C{1}{i} = C{2}{i}
        fid = fopen(fullfile(path, metaName), 'r');
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


% =========================================================
% Read nSamp timepoints from the binary file, starting
% at timepoint offset samp0. The returned array has
% dimensions [nChan,nSamp]. Note that nSamp returned
% is the lesser of: {nSamp, timepoints available}.
%
% IMPORTANT: samp0 and nSamp must be integers.
%
    function dataArray = ReadBin(samp0, nSamp, meta, binName, path)

        nChan = str2double(meta.nSavedChans);

        nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
        samp0 = max(samp0, 0);
        nSamp = min(nSamp, nFileSamp - samp0);

        sizeA = [nChan, nSamp];

        fid = fopen(fullfile(path, binName), 'rb');
        fseek(fid, samp0 * 2 * nChan, 'bof');
        dataArray = fread(fid, sizeA, 'int16=>double');
        fclose(fid);
    end % ReadBin


% =========================================================
% Return sample rate as double.
%
    function srate = SampRate(meta)
        if strcmp(meta.typeThis, 'imec')
            srate = str2double(meta.imSampRate);
        else
            srate = str2double(meta.niSampRate);
        end
    end % SampRate


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
% Having acquired a block of raw imec data using ReadBin(),
% convert values to gain-corrected voltages. The conversion
% is only applied to the saved-channel indices in chanList.
% Remember saved-channel indices are in range [1:nSavedChans].
% The dimensions of the dataArray remain unchanged. ChanList
% examples:
%
%   [1:AP]      % all AP chans (AP from ChannelCountsIM)
%   [2,6,20]    % just these three channels
%
    function dataArray = GainCorrectIM(dataArray, chanList, meta)

        % Look up gain with acquired channel ID
        chans = OriginalChans(meta);
        [APgain,LFgain] = ChanGainsIM(meta);
        nAP = length(APgain);
        nNu = nAP * 2;

        % Common conversion factor
        fI2V = Int2Volts(meta);

        for i = 1:length(chanList)
            j = chanList(i);    % index into timepoint
            k = chans(j);       % acquisition index
            if k <= nAP
                conv = fI2V / APgain(k);
            elseif k <= nNu
                conv = fI2V / LFgain(k - nAP);
            else
                continue;
            end
            dataArray(i,:) = dataArray(i,:) * conv;
        end
    end

    
end % ReadBin