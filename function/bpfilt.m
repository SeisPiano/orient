function filt_data = bpfilt(data,freqmin,freqmax,dt,corners)
% Bandpass filtering of time series.
% usage: filt_data = bpfilt(data,low frequency (Hz),high frequency (Hz),sampling interval (s))
%
% Yuechu Wu
% 12131066@mail.sustech.edu.cn
% 2022-10-08

if nargin == 4
    corners = 2; 
end

samprate = 1/dt;		   % sampling frequency

ch1n = 2*freqmin/samprate;
ch2n = 2*freqmax/samprate;

[b,a] = butter(corners,[ch1n ch2n],'bandpass');
filt_data = filtfilt(b,a,data);

return
