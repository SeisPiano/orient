function snr = calsnr(data,otime,btime,dist,dt,snr_op)
% original version was written by Dr. Ba Manh Le.
% the snr is calculated by sqrt(∑signal^2/∑noise^2)

if nargin == 5
    snr_op = 1;
end

if snr_op == 1

    rayl_vel = [3 4]; % lower and upper limits of Rayleigh wave velocity
    nois_vel = [2 2.5]; % lower and upper limits of noise velocity
    sigwin = round((dist./fliplr(rayl_vel)+otime-btime)/dt)+1;
    noiwin = round((dist./fliplr(nois_vel)+otime-btime)/dt)+1;

    if noiwin(2) > length(data)
        noiwin(2) = length(data);
    end
    if noiwin(1) >= length(data)
        snr = 0;
    else

        snr = sum(data(sigwin(1):sigwin(2)).^2)/sum(data(noiwin(1):noiwin(2)).^2);
        snr = sqrt(snr);
        
    end

elseif snr_op == 2

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

else
    error('Incorrect snr_op! Expected value is 1 or 2.');
end

return