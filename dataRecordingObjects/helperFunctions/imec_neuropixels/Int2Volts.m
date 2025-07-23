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
switch meta.typeThis
    case 'imec'
        if isfield(meta,'imMaxInt')
            maxInt = str2num(meta.imMaxInt);
        else
            maxInt = 512;
        end
        fI2V = str2double(meta.imAiRangeMax) / maxInt;
    case 'nidq'
        fI2V = str2double(meta.niAiRangeMax) / str2double(meta.niMaxInt);
    case 'obx'
        fI2V = str2double(meta.obAiRangeMax) / str2double(meta.obMaxInt);
end
end % Int2Volts

