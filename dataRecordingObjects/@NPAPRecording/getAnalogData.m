function [V_uV,t_ms]=getAnalogData(obj,channels,startTime_ms,window_ms)
%Extract Neuralynx recording data from file to memory
%Usage: [V_uV,t_ms]=obj.getData(channels,startTime_ms,window_ms);
%Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
%        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
%        window_ms - a scalar [1x1] with the window duration [ms].
%Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
%        t_ms - A time vector relative to recording start (t=0 at start)



%Implementation:
binName = string(obj.dataFileNames(contains(obj.dataFileNames,'nidq')));

path = obj.recordingDir;

meta = ReadMeta(binName, path);

chanList = channels;

nSamp = round((window_ms/1000)*str2double(meta.niSampRate));

nChan = str2double(meta.nSavedChans);

nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);

nSamp = min(nSamp, nFileSamp);


V_uV = zeros(length(channels),length(startTime_ms),nSamp);

origianlTimeStamps=1:0.05:5000;
onsetSync = interp1(originSQW'*1000, Neur'*1000, origianlTimeStamps, 'linear');


for trials = 1:length(startTime_ms)

    t0= startTime_ms(trials)/1000;

    samp0 = round(t0*str2double(meta.niSampRate));

    dataArray = ReadBin(samp0, nSamp, meta, binName, path, obj.convertData2Double);

    dataArray = dataArray(chanList,:);

    if obj.convertData2Double
        V_uV(:,trials,:) = GainCorrectNI(dataArray, chanList, meta); %convert to microvolts
    else
        V_uV = dataArray;
    end

end


%digits(6)
t_ms = (1:numel(V_uV))/obj.samplingFrequencyNI*1000;
%Added by Mark 28/6/25 - the code below resulted in size differences between V_uV and t_ms
%t_ms = 0:str2double(meta.niSampRate)/1000:str2double(meta.fileTimeSecs)*1000;
%t_ms= round(t_ms*1e6)/1e6;

end


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
function dataArray = ReadBin(samp0, nSamp, meta, binName, path, convertData2Double)

nChan = str2double(meta.nSavedChans);

nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
samp0 = max(samp0, 0);
nSamp = min(nSamp, nFileSamp - samp0);

sizeA = [nChan, nSamp];

fid = fopen(fullfile(path, binName), 'rb');
fseek(fid, samp0 * 2 * nChan, 'bof');

if convertData2Double
    dataArray = fread(fid, sizeA, 'int16=>double');
else
    dataArray = fread(fid, sizeA, 'int16');
end

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
% Having acquired a block of raw nidq data using ReadBin(),
% convert values to gain-corrected voltages. The conversion
% is only applied to the saved-channel indices in chanList.
% Remember saved-channel indices are in range [1:nSavedChans].
% The dimensions of the dataArray remain unchanged. ChanList
% examples:
%
%   [1:MN]      % all MN chans (MN from ChannelCountsNI)
%   [2,6,20]    % just these three channels
%
function dataArray = GainCorrectNI(dataArray, chanList, meta)

[MN,MA] = ChannelCountsNI(meta);
fI2V = Int2Volts(meta);

for i = 1:length(chanList)
    j = chanList(i);    % index into timepoint
    conv = fI2V / ChanGainNI(j, MN, MA, meta);
    dataArray(i,:) = dataArray(i,:) * conv;
end
end