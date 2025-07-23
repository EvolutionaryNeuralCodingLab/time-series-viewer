
% =========================================================
% Return sample rate as double.
%
function srate = SampRate(meta)
    switch meta.typeThis
        case 'imec'
            srate = str2double(meta.imSampRate);
        case 'nidq'
            srate = str2double(meta.niSampRate);
        case 'obx'
            srate = str2double(meta.obSampRate);
    end
end % SampRate