function filt_data = bpfilt(data,freqmin,freqmax,dt,corners,zerophase)
% Bandpass filtering of time series.
% usage: filt_data = bpfilt(data, low corner frequency (Hz), high corner frequency (Hz),
% sampling interval (s), filter corners / order, zerophase (true or false))
%
% Yuechu Wu
% 12131066@mail.sustech.edu.cn
% 2022-10-08

if nargin == 4
    corners = 2;
    zerophase = true;
end

if nargin == 5
    zerophase = true;
end

samprate = 1/dt;		   % sampling frequency

ch1n = 2*freqmin/samprate;
ch2n = 2*freqmax/samprate;

[b,a] = butter(corners,[ch1n ch2n],'bandpass');

if zerophase
    filt_data = filtfilt(b,a,data);
else
    filt_data = filter(b,a,data);
end

return
