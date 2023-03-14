% plot_seismogram
%
% Plot three-component seismic waveform with surface wave of Z component 
% signal-to-noise ratio greater than threshold.
%
% The first 14 digits of the event file name should preferably be 
% yyyymmddHHMMSS.
%
% Yuechu Wu
% 12131066@mail.sustech.edu.cn
% 2023-03-09



clear; close all;

addpath ('function'); % path of matlab functions

network = 'XF';     % network name
stations = {'B01'}; % station names
ext = 'SAC';        % name of the event file extension

bhz = 'HHZ'; % Z component channel name
bh1 = 'HH1'; % 1 component channel name
bh2 = 'HH2'; % 2 component channel name

snr_op = 1;  % options for using SNR calculations: (1) sqrt(∑signal^2/∑noise^2)
%                                                  (2) max(abs(signal)/mean(abs(noise))
minsnr = 5;  % signal to noise ratio threshold

filtband = [0.02 0.05]; % bandpass filter band in Hz
tap_width = 0.1;        % width in percent of attenuation window

INPUTdir  = 'DATA/EVENT';  % directory for data input
FIGUREdir = 'FIGURES';     % directory for figure output



%%%%% END OF USER INPUT %%%%%
channels = {bhz,bh1,bh2};
SELECTdir = sprintf('%s/SEISMOGRAMS/SELECT',FIGUREdir);

if ~exist(FIGUREdir,'dir')
    mkdir(FIGUREdir);
end
if ~exist(SELECTdir,'dir')
    mkdir(SELECTdir);
end
if ~exist(fullfile(SELECTdir,network),'dir')
    mkdir(fullfile(SELECTdir,network));
end

filesuff = sprintf('*%s*%s',bhz,ext); % match the Z component event file name

for ista = 1:length(stations) % begin station loop
    station = stations{ista};

    if ~exist(fullfile(SELECTdir,network,station),'dir')
        mkdir(fullfile(SELECTdir,network,station));
    end

    figoutpath = sprintf('%s/%s/%s',SELECTdir,network,station);

    inpath = sprintf('%s/%s/%s',INPUTdir,network,station);
    filenames = dir(fullfile(inpath,filesuff));

    for ie = 1:length(filenames)  % begin event loop

        close all; clear traces

        filename_z = sprintf('%s/%s',inpath,filenames(ie).name);

        eventid = filenames(ie).name(1:14); % yyyymmddHHMMSS
        event_otime = sprintf('%s-%s-%s %s:%s',eventid(1:4),eventid(5:6),eventid(7:8),eventid(9:10),eventid(11:12));
        disp(event_otime); % yyyy-mm-dd HH:MM

        filename_1 = strrep(filename_z,bhz,bh1);
        filename_2 = strrep(filename_z,bhz,bh2);

        if ~exist(filename_1,'file') || ~exist(filename_2,'file')
            disp('Skipping. At least one horizontal component is missing');
            continue
        end

        rawdata_z = rsac(filename_z);

        [otime_z,btime_z,dist,baz,dt] = lh(rawdata_z,'O','B','DIST','BAZ','DELTA');
        [gcarc,magnitude] = lh(rawdata_z,'GCARC','MAG');
        data_z = datapreproc(rawdata_z(:,2),dt,'filtband',filtband,'tap_width',tap_width);

        snr_z = calsnr(data_z,otime_z,btime_z,dist,dt,snr_op);

        if snr_z < minsnr
            disp('Skipping. Low SNR');
            continue
        end

        rawdata_1 = rsac(filename_1);
        rawdata_2 = rsac(filename_2);

        data_length = [length(rawdata_z(:,2)),length(rawdata_1(:,2)),length(rawdata_2(:,2))];

        if length(unique(data_length))>1
            disp('Skipping. Different data lengths')
            continue
        end

        data_1 = datapreproc(rawdata_1(:,2),dt,'filtband',filtband,'tap_width',tap_width);
        data_2 = datapreproc(rawdata_2(:,2),dt,'filtband',filtband,'tap_width',tap_width);

        event_length = length(data_z)*dt;
        time = transpose(0:dt:event_length-dt);

        traces(1).data = data_z;
        traces(2).data = data_1;
        traces(3).data = data_2;

        figure(1)
        set(gcf,'PaperPositionMode','manual');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperOrientation','portrait');
        set(gcf,'PaperPosition',[.05 .05 8 15]);

        for ic = 1:length(channels)
            subplot(length(channels),1,ic);
            plot(time,traces(ic).data,'color','#000000','LineWidth',0.5);
            title(channels{ic})
            xlim([0 event_length])
            ylabel('Displacement (m)')
            xlabel('Time (s)')
        end

        event_info = sprintf('%s, Distance=%0.1f°, M=%0.1f',event_otime,gcarc,magnitude);
        sgtitle(sprintf('%s %s',station,event_info));

        filename = sprintf('%s/%s_%s_%s.png',figoutpath,eventid,network,station);
        print('-dpng',filename);


    end % end event loop

end % end station loop
