function snr = calsnr2(data,otime,btime,dist,dt)
% This function refers to the noise cross-correlation method of Huajian Yao
% Reference: http://yaolab.ustc.edu.cn/resources.php?i=28
% the snr is calculated by max(abs(signal)/mean(abs(noise))
%
% 2022-10-25
% Yuechu Wu
% 12131066@mail.sustech.edu.cn

velovity = [0.5, 5]; % lower and upper limits of Rayleigh wave velocity
noise_win_len = 150; % length of noise window (s), which is defined right after the signal window

signal_win = round((dist./fliplr(velovity)+otime-btime)/dt)+1;
noise_win  = [signal_win(2), signal_win(2) + round(noise_win_len/dt)];

if signal_win(1) < 1
    signal_win(1) = 1;
end

if noise_win(2) > length(data)
    noise_win(2) = length(data);
end

if noise_win(1) >= length(data)
    snr = 0;
else
    data_enve = abs(hilbert(data));
    signal_amp_max = max(data_enve(signal_win));
    noise_amp_mean = mean(data_enve(noise_win));

    if noise_amp_mean == 0
        if signal_amp_max > 0
            snr = 100;
        else
            snr = 0;
        end
    else
        snr = signal_amp_max/noise_amp_mean; % max signal envelope amplitue / average noise envelope amplitude
    end

end

return