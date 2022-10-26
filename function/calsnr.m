function snr = calsnr(data,otime,btime,dist,dt)
% original version was written by Dr. Ba Manh Le.
% the snr is calculated by sqrt(∑signal^2/∑noise^2)
rayl_vel = [3 4]; % lower and upper limits of Rayleigh wave velocity
nois_vel = [2 2.5]; % lower and upper limits of noise velocity
sigwin = round((dist./fliplr(rayl_vel)+otime-btime)/dt)+1;
noiwin = round((dist./fliplr(nois_vel)+otime-btime)/dt)+1;
if sigwin(1)<1
    sigwin(1)=1;
end
if noiwin(2)>length(data)
    noiwin(2)=length(data);
end
if noiwin(1)>=length(data) || sigwin(2)<1
    snr = 0;
else
    lensigwin = sigwin(2)-sigwin(1);
    lennoiwin = noiwin(2)-noiwin(1);

    if lensigwin > lennoiwin
        snr = 0;
    else
        snr = sum(data(sigwin(1):sigwin(2)).^2)/sum(data(noiwin(1):noiwin(2)).^2);
        snr = sqrt(snr);
    end
end

return