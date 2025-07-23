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

    np1_imro = [0,1020,1030,1200,1100,1120,1121,1122,1123,1300];
    
    if isfield(meta,'imDatPrb_type')
        probeType = str2double(meta.imDatPrb_type);
    else
        probeType = 0;
    end

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
        elseif k <= nNu && ismember(probeType, np1_imro) 
            %If it is NP1.0 and k is not sync channel, calculate conversion
            %for LFP channels 
            conv = fI2V / LFgain(k - nAP);
        else
            continue;
        end
        dataArray(j,:) = dataArray(j,:) * conv;
    end
end % GainCorrectIM