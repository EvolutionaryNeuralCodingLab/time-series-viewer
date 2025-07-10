
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
