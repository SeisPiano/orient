function filt_data = bpfilt(data,dt,lowfreq,highfreq)
% Bandpass filtering of time series.
% usage: filt_data = bpfilt(data, sampling interval [s], low frequency [Hz], high frequency [Hz])
%
% Yuechu Wu
% 12131066@mail.sustech.edu.cn
% 2022-10-08

fn = 1/2/dt;
[b,a] = butter(2,[lowfreq/fn,highfreq/fn]);
filt_data = filtfilt(b,a,data);

return
