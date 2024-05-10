% plot_waveform
%
% Plot three-component seismic waveform with surface wave of Z component 
% signal-to-noise ratio greater than threshold.
%
% The first 15 digits of the event file name should preferably be 
% yyyymmdd_HHMMSS.
%
% Yuechu Wu
% 12131066@mail.sustech.edu.cn
% 2023-03-09
%
% Update to adapt to the new version of DownloadSeisData
% Yuechu Wu


clear; close all;

addpath ('function');  % path of matlab functions

INPUTdir  = 'DATA/sacdata_event';  % directory for data input
FIGUREdir = 'FIGURES';             % directory for figure output

network = 'XF';      % network name
stations = {'B01'};  % station names
ext = 'SAC';         % name of the event file extension

channel_z = 'HHZ';  % Z component channel name
channel_1 = 'HH1';  % 1 component channel name
channle_2 = 'HH2';  % 2 component channel name

snr_op = 1;     % options for using SNR calculations. 1: sqrt(∑signal^2/∑noise^2)
                %                                     2: max(abs(signal)/mean(abs(noise))

minsnr = 8;            % signal to noise ratio threshold

filt_band = [0.03, 0.04];  % bandpass filter band in Hz
tap_width = 0.1;           % width in percent of attenuation window
velocity = [3, 4.5];       % lower and upper limits of Rayleigh wave velocity





%%%%% END OF USER INPUT %%%%%
channels = {channel_z,channel_1,channle_2};
SELECTdir = sprintf('%s/waveform/select',FIGUREdir);

if ~exist(FIGUREdir,'dir')
    mkdir(FIGUREdir);
end

if ~exist(SELECTdir,'dir')
    mkdir(SELECTdir);
end

if ~exist(fullfile(SELECTdir,network),'dir')
    mkdir(fullfile(SELECTdir,network));
end


% List the contents of the folder and remove the . and ..
eventids = setdiff({dir(INPUTdir).name},{'.','..'});


for ista = 1:length(stations) % begin station loop
    station = stations{ista};
    disp(station)
    close all; clear st

    if ~exist(fullfile(SELECTdir,network,station),'dir')
        mkdir(fullfile(SELECTdir,network,station));
    end

    figoutpath = sprintf('%s/%s/%s',SELECTdir,network,station);


    for ie = 1:length(eventids) % begin event loop
        close all
        eventid = eventids{ie};
        event_otime = sprintf('%s-%s-%sT%s:%s:%s',eventid(1:4),eventid(5:6),eventid(7:8), ...
                      eventid(10:11),eventid(12:13),eventid(14:15));
        disp(event_otime); % yyyy-mm-dd HH:MM:SS


        % match the Z component event file name
        file_z = dir(sprintf('%s/%s/*%s*%s*%s',INPUTdir,eventid,station,channel_z,ext));
        if isempty(file_z)
            continue
        end
        
        filename_z = sprintf('%s/%s/%s',INPUTdir,eventid,file_z.name);

        filename_1 = strrep(filename_z,channel_z,channel_1);
        filename_2 = strrep(filename_z,channel_z,channle_2);

        if ~exist(filename_1,'file') || ~exist(filename_2,'file')
            disp('Skipping. At least one horizontal component is missing');
            continue
        end

        tr_z = rsac(filename_z);

        [otime_z,btime_z,dist,dt] = lh(tr_z,'O','B','DIST','DELTA');
        [gcarc,magnitude] = lh(tr_z,'GCARC','MAG');
        data_z = datapreproc(tr_z(:,2),dt,'filt_band',filt_band,'tap_width',tap_width);

        snr_z = calsnr(data_z,otime_z,btime_z,dist,dt,snr_op);

        if snr_z < minsnr
            disp('Skipping. Low SNR');
            continue
        end

        tr_1 = rsac(filename_1);
        tr_2 = rsac(filename_2);

        % data_length = [length(tr_z(:,2)),length(tr_n(:,2)),length(tr_e(:,2))];

        % if length(unique(data_length)) > 1
        %     disp('Skipping. Different data lengths')
        %     continue
        % end

        data_1 = datapreproc(tr_1(:,2),dt,'filt_band',filt_band,'tap_width',tap_width);
        data_2 = datapreproc(tr_2(:,2),dt,'filt_band',filt_band,'tap_width',tap_width);

        event_length = length(data_z)*dt;
        time = transpose(btime_z:dt:event_length+btime_z-dt);
        signal_window = dist./fliplr(velocity);

        st(1).data = data_z;
        st(2).data = data_1;
        st(3).data = data_2;

        figure(1)
        set(gcf,'PaperPositionMode','manual');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperOrientation','portrait');
        set(gcf,'PaperPosition',[.05 .05 8 15]);

        for ic = 1:length(channels)
            subplot(length(channels),1,ic);
            plot(time,st(ic).data,'color','#000000','LineWidth',0.5);hold on
            xline(signal_window(1),'color','#FF8C00');
            xline(signal_window(2),'color','#FF8C00');
            title(channels{ic})
            xlim([btime_z event_length+btime_z])
            ylabel('Displacement (m)')
            xlabel('Time (s)')
        end

        event_info = sprintf('%s, Distance=%0.1f°, M=%0.1f',event_otime,gcarc,magnitude);
        sgtitle(sprintf('%s %s',station,event_info));

        filename = sprintf('%s/%s_%s_%s.png',figoutpath,eventid,network,station);
        print('-dpng',filename);


    end % end event loop

end % end station loop
