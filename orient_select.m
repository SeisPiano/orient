% orient_select
%
% Determination of ocean bottom seismograph (OBS) orientation via
% Rayleigh-wave polarization.
%
% The script orient_select has the same function of the script orient, 
% but it reduce the error by using selected events.
% The selected event seismograms should be manually selected after 
% running plot_seismogram, which in the SELECTdir path.
%
% Yuechu Wu
% 2023-03-12
% 12131066@mail.sustech.edu.cn
%
% Update to adapt to the new version of DownloadSeisData
% Yuechu Wu

clear; close all;

addpath ('function');  % path of matlab functions

network = 'XF';      % network name
stations = {'B01'};  % station names
ext = 'SAC';         % name of the event file extension

% If you don't know which horizontal channel corresponds to the N component,
% try it. Only the correct channels can get reasonable results.
channel_z = 'HHZ';  %  Z  component channel name
channel_n = 'HH2';  % nominal N component channel name
channel_e = 'HH1';  % nominal E component channel name

bpass = [0.015, 0.025, 0.035, 0.045];  % frequency used to calculate
dbp = 0.005;           % the band for each frequency bpass +/- dbp
tap_width = 0.1;       % width in percent of attenuation window
velocity = [3, 4.5];   % lower and upper limits of Rayleigh wave velocity

mincorr = 0.8;  % minimum cross correlation use to pick the result
snr_op = 1;     % options for using SNR calculations. 1: sqrt(∑signal^2/∑noise^2)
                %                                     2: max(abs(signal)/mean(abs(noise))

minsnr = 8;            % signal to noise ratio threshold
min_similarity = 0.8;  % minimum similar amplitude magnitude
isoverwrite = 1;       % overwrite output files

INPUTdir = 'DATA/sacdata_event';  % directory for event data
OUTPUTdir = 'DATA/ORIENT';        % directory for data output
FIGUREdir = 'FIGURES';            % directory for figure output


%%%%% END OF USER INPUT %%%%%
SELECTdir = sprintf('%s/SEISMOGRAMS/SELECT',FIGUREdir);

if ~exist(OUTPUTdir,'dir')
    mkdir(OUTPUTdir);
end

if ~exist(FIGUREdir,'dir')
    mkdir(FIGUREdir);
end

if ~exist(fullfile(OUTPUTdir,network),'dir')
    mkdir(fullfile(OUTPUTdir,network));
end

if ~exist(fullfile(FIGUREdir,network),'dir')
    mkdir(fullfile(FIGUREdir,network));
end

figoutpath = sprintf('%s/%s',FIGUREdir,network);


for ista = 1:length(stations) % begin station loop
    station = stations{ista};
    disp(station);

    outfile = sprintf('%s/%s/%s_%s_orient_select.txt',OUTPUTdir,network,network,station);

    if ~exist(outfile,'file') || isoverwrite

        close all; clear outdata iee;

        iee = 0; % effective event

        selectpath = sprintf('%s/%s/%s',SELECTdir,network,station);
        select_filenames = dir(fullfile(selectpath,'*.png'));

        for ie = 1:length(select_filenames) % begin event loop

            eventid = select_filenames(ie).name(1:15); % yyyymmdd_HHMMSS
            event_otime = sprintf('%s-%s-%sT%s:%s:%s',eventid(1:4),eventid(5:6),eventid(7:8), ...
                          eventid(10:11),eventid(12:13),eventid(14:15));
            disp(event_otime); % yyyy-mm-dd HH:MM:SS

            
            % match the Z component event file name
            file_z = dir(sprintf('%s/%s/*%s*%s*%s',INPUTdir,eventid,station,channel_z,ext));

            if isempty(file_z)
                continue
            end

            filename_z = sprintf('%s/%s/%s',INPUTdir,eventid,file_z.name);

            filename_n = strrep(filename_z,channel_z,channel_n);
            filename_e = strrep(filename_z,channel_z,channel_e);

            if ~exist(filename_n,'file') || ~exist(filename_e,'file')
                disp('Skipping. At least one horizontal component is missing');
                continue
            end

            tr_z = rsac(filename_z);
            tr_n = rsac(filename_n);
            tr_e = rsac(filename_e);

            % data_length = [length(tr_z(:,2)),length(tr_n(:,2)),length(tr_e(:,2))];

            % if length(unique(data_length)) > 1
            %     disp('Skipping. Different data lengths')
            %     continue
            % end

            [otime_z,btime_z,dist,baz,dt] = lh(tr_z,'O','B','DIST','BAZ','DELTA');
            [otime_n,btime_n] = lh(tr_n,'O','B');
            [otime_e,btime_e] = lh(tr_e,'O','B');

            [gcarc,magnitude] = lh(tr_z,'GCARC','MAG');

            for ifreq = 1:length(bpass) % begin frequency loop

                filt_band = [bpass(ifreq)-dbp, bpass(ifreq)+dbp];
                fprintf('%.0f-%.0f s\n',1/(bpass(ifreq)+dbp),1/(bpass(ifreq)-dbp));
           
                data_z = datapreproc(tr_z(:,2),dt,'filt_band',filt_band,'tap_width',tap_width);

                snr_z = calsnr(data_z,otime_z,btime_z,dist,dt,snr_op);

                if snr_z < minsnr
                    disp('Skipping. Low SNR');
                    continue
                end
     
                data_n = datapreproc(tr_n(:,2),dt,'filt_band',filt_band,'tap_width',tap_width);       
                data_e = datapreproc(tr_e(:,2),dt,'filt_band',filt_band,'tap_width',tap_width);

                sigwin_z = round((dist./fliplr(velocity)+otime_z-btime_z)/dt)+1;
                sigwin_n = round((dist./fliplr(velocity)+otime_n-btime_n)/dt)+1;
                sigwin_e = round((dist./fliplr(velocity)+otime_e-btime_e)/dt)+1;

                data_z_hilb = imag(hilbert(data_z));

                if length(data_e)>=sigwin_e(2) && length(data_n)>=sigwin_n(2) && length(data_z_hilb)>=sigwin_z(2)
                    sig_z = data_z_hilb(sigwin_z(1):sigwin_z(2));
                    sig_n = data_n(sigwin_n(1):sigwin_n(2));
                    sig_e = data_e(sigwin_e(1):sigwin_e(2));
                else
                    continue
                end

                % only use events with similar amplitude magnitude
                if sum(sig_z.^2)/sum(sig_n.^2) >= min_similarity && sum(sig_z.^2)/sum(sig_e.^2) >= min_similarity
                    sig_z = sig_z/max(sig_z);
                    sig_n = sig_n/max(sig_n);
                    sig_e = sig_e/max(sig_e);
                else
                    continue
                end

                if length(sig_z) == length(sig_n) && length(sig_z) == length(sig_e)

                    [corr,x,corr2,x2] = calcOrient(sig_z,sig_n,sig_e,baz);

                    iee = iee + 1; % effective event + 1

                    % each column of outdata represents the effective earthquake serial number
                    % Each row of outdata represents the identifier of frequency, identifier of
                    % event, back azimuth, cross correlation, x, second cross correlation and x2.
                    outdata(iee,:) = [ifreq, ie, baz, corr, x, corr2, x2, snr_z];
                else
                    continue
                end


            end % end frequency loop

        end % end event loop

        save(outfile,'outdata','-ascii');

    end % end overwrite


    plotOrient2(outfile,mincorr,network,station); % network and station are not necessary

    figure(1)
    set(gcf,'PaperPositionMode','manual');
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperOrientation','portrait');
    set(gcf,'PaperPosition',[.05 1.5 8 4]);
    xlabel('Correlation');
    ylabel('Degrees from N'); % Degrees from the direction of N component
    pdffile1 = sprintf('%s/%s_%s_north_select.pdf',figoutpath,network,station);
    print('-dpdf',pdffile1);

    figure(2)
    pdffile2 = sprintf('%s/%s_%s_rose_select.pdf',figoutpath,network,station);
    print('-dpdf',pdffile2);

end % end station loop
