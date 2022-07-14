function writeKwik( filename, data, triggers, triggersChNum, ...
  recordingDuration, samplingFreq, bitVolts )
%WRITEKWIK Writes two kwik files containing data and triggers
%   This function creates two files following the KWIK specifications
%   (here: https://github.com/klusta-team/kwiklib/wiki/Kwik-format#kwik)
%   
%   The file filename.kwd contains the raw data, while the file
%   filename.kwe is the file containing the trigger information. Two files
%   are used following the variations of the original format specified
%   here: https://open-ephys.atlassian.net/wiki/display/OEW/KWIK+format
%
%   INPUT: 
%     filename the desired filename 
%     data the raw data, passed as a 3D matrix of the format [nChannels x nTrials x nSamples]
%     triggers the triggers information, passed as cell array
%     triggersChNum channel number information for the triggers
%     recordingDuration duration of recordings
%     samplingFreq
%     bitVolts microvoltsPerAD
%   OUTPUT: none. 2 files will be written on disk
%
%   Author: Stefano.Masneri@brain.mpg.de
%   Date: 18.11.2016

dataFile = [filename '.kwd'];
trigFile = [filename '.kwe'];
numCh = size(data, 1);
numSamples = size(data, 3);

if isscalar(samplingFreq)
  samplingFreq = single(repmat(samplingFreq, [1, numCh]));
end

if isscalar(bitVolts)
  bitVolts = single(repmat(bitVolts, [1, numCh]));
end

datatype = class(data);
if strcmpi('uint16', datatype) || strcmpi('int16', datatype)
  bitDepth = 16;
elseif strcmpi('uint8', datatype) || strcmpi('int8', datatype)
  bitDepth = 8;
elseif strcmpi('uint32', datatype) || strcmpi('int32', datatype)
  bitDepth = 32;
else
  bitDepth = 64;
end
bitDepth = int16(bitDepth);

% first write the data in the .kwd file
h5create(dataFile, '/recordings/0/data', size(squeeze(data(:, 1, :))), ...
  'Datatype', class(data));
h5write(dataFile, '/recordings/0/data', squeeze(data(:,1,:)));

% then write the info in application_data
h5create(dataFile, '/recordings/0/application_data/channel_bit_volts', size(bitVolts));
h5write(dataFile, '/recordings/0/application_data/channel_bit_volts', bitVolts);

h5create(dataFile, '/recordings/0/application_data/channel_sample_rates', size(samplingFreq));
h5write(dataFile, '/recordings/0/application_data/channel_sample_rates', samplingFreq);

% write attribute of recording
h5writeatt(dataFile, '/recordings/0', 'sample_rate', samplingFreq(1)); 
h5writeatt(dataFile, '/recordings/0', 'bit_depth', bitDepth);
h5writeatt(dataFile, '/recordings/0', 'name', filename); 
h5writeatt(dataFile, '/recordings/0', 'start_sample', 0);
h5writeatt(dataFile, '/recordings/0', 'start_time', 0);

% write attribute of application_data
isMultiSampleRateData = 0;
if length(unique(samplingFreq)) > 1
  isMultiSampleRateData = 1;
end
h5writeatt(dataFile, '/recordings/0/application_data', 'is_MultiSampleRate_data', isMultiSampleRateData); 

%write attribute of data
validSamples = int32(repmat(numSamples, [numCh, 1]));
h5writeatt(dataFile, '/recordings/0/data', 'valid_samples', validSamples);

% for k = 1:numCh
%   h5create(dataFile, ['/recordings/' num2str(k) '/data'], size(data(k, :, :)), ...
%     'Datatype', class(data));
%   h5write( dataFile, ['/recordings/' num2str(k) '/data'], data(k, :, :));
%   h5writeatt(dataFile, ['/recordings/' num2str(k)], 'sample_rate', samplingFreq(k)); 
%   h5writeatt(dataFile, ['/recordings/' num2str(k)], 'bit_depth', bitDepth);
%   % now write the attribute in application data
%   % require the low level functions
%   plist = 'H5P_DEFAULT';
%   fid = H5F.open(dataFile, 'H5F_ACC_RDWR', plist);
%   gid = H5G.create(fid, ['/recordings/' num2str(k) '/application_data'], plist, plist, plist);
%   H5G.close(gid);
%   H5F.close(fid);
%   h5writeatt(dataFile, ['/recordings/' num2str(k) '/application_data'], 'channel_bit_volts', bitVolts(k));
% end

% write some global attributes
h5writeatt(dataFile, '/', 'datatype', datatype);
h5writeatt(dataFile, '/', 'kwik_version', 2);
h5writeatt(dataFile, '/', 'recordingDuration', recordingDuration);

% now for the events
lengthsTriggers = cellfun('length', triggers); % get the length of each element of the cell
trigMat = cell2mat(triggers);                  % concatenate all the cell elements
h5create(trigFile, '/event_types/TTL/events/time_samples', size(trigMat), ...
  'Datatype', class(trigMat));
h5write(trigFile, '/event_types/TTL/events/time_samples', trigMat);

% set the on off value for the triggers
triggersOnOff = zeros(1, length(trigMat));     
triggersOnOff(1:lengthsTriggers(1:2:end)) = 1;
triggersOnOff = uint8(triggersOnOff);
h5create(trigFile, '/event_types/TTL/events/user_data/eventID', size(triggersOnOff), ...
  'Datatype', class(triggersOnOff));
h5write(trigFile, '/event_types/TTL/events/user_data/eventID', triggersOnOff);

% set the channel value for the triggers
channelInfo = [];
for k = 1:length(triggersChNum)
  channelInfo(end+1 : end+lengthsTriggers(k)) = triggersChNum(k);
end
h5create(trigFile, '/event_types/TTL/events/user_data/event_channels', size(channelInfo), ...
  'Datatype', class(channelInfo));
h5write(trigFile, '/event_types/TTL/events/user_data/event_channels', channelInfo);

end