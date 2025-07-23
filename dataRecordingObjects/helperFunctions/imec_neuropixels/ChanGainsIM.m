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

% =========================================================
% Return gain arrays for imec channels.
%
% Index into these with original (acquired) channel IDs.
%
function [APgain,LFgain, APChan0_to_uV, LFChan0_to_uV] = ChanGainsIM(meta)
    % list of probe types with NP 1.0 imro format
    np1_imro = [0,1020,1030,1200,1100,1120,1121,1122,1123,1300];
    % number of channels acquired
    acqCountList = str2num(meta.acqApLfSy);

    APgain = zeros(acqCountList(1));     % default type = float64
    LFgain = zeros(acqCountList(2));     % empty array for 2.0

    if isfield(meta,'imDatPrb_type')
        probeType = str2double(meta.imDatPrb_type);
    else
        probeType = 0;
    end

    if ismember(probeType, np1_imro)
        % imro + probe allows setting gain independently for each channel
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
    else               
        % get gain from  imChan0apGain, if present
        if isfield(meta,'imChan0apGain')
            APgain = APgain + str2num(meta.imChan0apGain);
            if acqCountList(2) > 0               
                LFgain = LFgain + str2num(meta.imChan0lfGain);
            end
        elseif (probeType == 1110)
            % active UHD, for metadata lacking imChan0apGain, get gain from
            % imro table header
            currList = sscanf(meta.imroTbl, '(%d,%d,%d,%d,%d');
            APgain = APgain + currList(4);
            LFgain = LFgain + currList(5);
        elseif (probeType == 21) || (probeType == 24)
            % development NP 2.0; APGain = 80 for all AP
            % return 0 for LFgain (no LF channels)
            APgain = APgain + 80;        
        elseif (probeType == 2013)
            % commercial NP 2.0; APGain = 80 for all AP
            APgain = APgain + 100;
        else
            fprintf('unknown gain, setting APgain to 1\n');
            APgain = APgain + 1;
        end
    end
    fI2V = Int2Volts(meta);
    APChan0_to_uV = 1e6*fI2V/APgain(1);
    if size(LFgain) > 0
        LFChan0_to_uV = 1e6*fI2V/LFgain(1);
    else
        LFChan0_to_uV = 0;
    end
end % ChanGainsIM