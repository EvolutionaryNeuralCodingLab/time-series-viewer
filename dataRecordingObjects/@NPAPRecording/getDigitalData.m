function [D,T_ms]=getDigitalData(obj,startTime_ms,window_ms)
%Extract MCRack digital data from file to memory
%Usage: [D,T_ms]=getDigitalData(startTime_ms,window_ms,name)
%Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
%        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
%        window_ms - a scalar [1x1] with the window duration [ms].
%        name - the name of the stream (if not entered, default name is used)
%Output: D - A 3D matrix [nChannels x nTrials x nSamples] with digitalData waveforms across specified channels and trials
%        T_ms - A time vector relative to recording start (t=0 at start)



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
% Return an array [lines X timepoints] of uint8 values for
% a specified set of digital lines.
%
% - dwReq is the one-based index into the saved file of the
%    16-bit word that contains the digital lines of interest.
% - dLineList is a zero-based list of one or more lines/bits
%    to scan from word dwReq.
%
    function digArray = ExtractDigital(dataArray, meta, dwReq, dLineList)
        % Get channel index of requested digital word dwReq
        if strcmp(meta.typeThis, 'imec')
            [AP, LF, SY] = ChannelCountsIM(meta);
            if SY == 0
                fprintf('No imec sync channel saved\n');
                digArray = [];
                return;
            else
                digCh = AP + LF + dwReq;
            end
        else
            [MN,MA,XA,DW] = ChannelCountsNI(meta);
            if dwReq > DW
                fprintf('Maximum digital word in file = %d\n', DW);
                digArray = [];
                return;
            else
                digCh = MN + MA + XA + dwReq;
            end
        end
        [~,nSamp] = size(dataArray);
        digArray = zeros(numel(dLineList), nSamp, 'uint8');
        for i = 1:numel(dLineList)
            digArray(i,:) = bitget(dataArray(digCh,:), dLineList(i)+1, 'int16');
        end
    end % ExtractDigital

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
    
    

    if startTime_ms == "Inf"
        nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
        


    end
    
    binName = string(obj.dataFileNames(contains(obj.dataFileNames,'nidq')));

    path = obj.recordingDir;

    meta = ReadMeta(binName, path);

    dw = 8;
    dLineList = [1, 2, 3];


end