% =========================================================
% Return gain for ith channel stored in the nidq file.
%
% ichan is a saved channel index, rather than an original
% (acquired) index.
%
% Note: there is no matching function for OBX, because
% the gain is fixed at 1.
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
