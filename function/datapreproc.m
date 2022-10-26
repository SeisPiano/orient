function preproc_data = datapreproc(raw_data,dt,filtfreq,tap_width)

% Data preprocessing
% usage: dataout = datapreproc(datain, sampling interval, filtering frequency band, taper width)
% tap_width can be 0, then the taper operation will not be executed.
%
% Yuechu Wu
% 12131066@mail.sustech.edu.cn
% 2022-09-29

rmean_data = raw_data - mean(raw_data);
rtrend_data = detrend(rmean_data);
filt_data = bpfilt(rtrend_data,dt,filtfreq(1),filtfreq(2));

if tap_width == 0
    preproc_data = filt_data;
else
    taper_data = filt_data.*tukeywin(length(filt_data),tap_width);
    preproc_data = taper_data;
end

return