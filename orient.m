% Determination of ocean bottom seismograph (OBS) orientation via
% Rayleigh-wave polarization.
%
% original version was written by Dr. Ba Manh Le.
% The input earthquake events must include the information of epicenter
% distance and back azimuth in SAC header variables.
%
% -- handle multiple stations;
% -- simplify multiple steps into functions;
% -- a method for calculating the signal-to-noise ratio has been added.
% Updated 2022-10-25
% Yuechu Wu
% 12131066@mail.sustech.edu.cn



clear; close all;

addpath ('function'); % path for matlab functions

network = 'XF'; % network name
stations = {'B02'}; % station names
ext = 'sac.disp';  % name of the extension

INPUTdir = 'DATA/EVENT'; % directory for data input
OUTPUTdir = 'DATA/ORIENT'; % directory for data output
FIGUREdir = 'FIGURES'; % directory for figure output

% If you don't know which horizontal channel corresponds to the N component,
% try it. Only the correct channels can get reasonable results.
bhz = 'HHZ'; %  Z  component channel name
bhn = 'HH2'; % N-S component channel name
bhe = 'HH1'; % E-W component channel name

bpass = [0.015, 0.025, 0.035, 0.045];  % frequency used to calculate
dbp = 0.005;  % the band for each frequency bpass +/- dbp
tap_width = 0.1; % width in percent of frequency vector of cosine taper
velocity = [3, 4.5]; % lower and upper limits of Rayleigh wave velocity

mincorr = 0.8;  % minimum cross correlation use to pick the result
snr_op = 2; % options for using SNR calculations: (1) sqrt(∑signal^2/∑noise^2)
%                                                 (2) max(abs(signal)/mean(abs(noise))
minsnr = 3;  % signal to noise ratio threshold
isoverwrite = 1; % overwrite output files

% END OF USER INPUT %
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

inpath = sprintf('%s/%s',INPUTdir,network);
outpath = sprintf('%s/%s',OUTPUTdir,network);
figoutpath = sprintf('%s/%s',FIGUREdir,network);

filesuff = sprintf('*%s*%s',bhz,ext); % change this for matching with file name format

for ista = 1:length(stations) % begin station loop
    station = stations{ista};
    outfile = sprintf('%s/%s_%s_orient.txt',outpath,network,station);

    if ~exist(outfile,'file') || isoverwrite

        close all; clear outdata iee;

        iee = 0; % effective event

        for ifreq = 1:length(bpass) % begin frequency loop
            lowfreq = bpass(ifreq)-dbp; higfreq = bpass(ifreq)+dbp;
            filenames = dir(sprintf('%s/%s/%s/%s',INPUTdir,network,station,filesuff));

            for ie = 1:length(filenames)  % begin event loop
                filename_z = sprintf('%s/%s/%s',inpath,station,filenames(ie).name);
                filename_n = strrep(filename_z,bhz,bhn);
                filename_e = strrep(filename_z,bhz,bhe);

                if exist(filename_n,'file') && exist(filename_e,'file')

                    fprintf('loading %s\n',filename_z);
                    rawdata_z = rsac(filename_z);

                    [otime_z,btime_z,dist,baz,dt] = lh(rawdata_z,'O','B','DIST','BAZ','DELTA');
                    data_z = datapreproc(rawdata_z(:,2),dt,[lowfreq,higfreq],tap_width);

                    snr_z = calsnr(data_z,otime_z,btime_z,dist,dt,snr_op);

                    dataz_hilb = imag(hilbert(data_z));

                    fprintf('loading %s\n',filename_n);
                    rawdata_n = rsac(filename_n);
                    [otime_n,btime_n,~,~,dt] = lh(rawdata_n,'O','B','DIST','BAZ','DELTA');
                    data_n = datapreproc(rawdata_n(:,2),dt,[lowfreq,higfreq],tap_width);
                    % snr_n = calsnr(datan,otimen,btimen,dist,dt,snr_op);

                    fprintf('loading %s\n',filename_e);
                    rawdata_e = rsac(filename_e);
                    [otime_e,btime_e,~,~,dt] = lh(rawdata_e,'O','B','DIST','BAZ','DELTA');
                    data_e = datapreproc(rawdata_e(:,2),dt,[lowfreq,higfreq],tap_width);
                    % snr_e = calsnr(datae,otimee,btimee,dist,dt,snr_op);


                    if snr_z > minsnr
                        sigwin_z = round((dist./fliplr(velocity)+otime_z-btime_z)/dt)+1;
                        sigwin_n = round((dist./fliplr(velocity)+otime_n-btime_n)/dt)+1;
                        sigwin_e = round((dist./fliplr(velocity)+otime_e-btime_e)/dt)+1;

                        if length(data_e)>=sigwin_e(2) && length(data_n)>=sigwin_n(2) && length(dataz_hilb)>=sigwin_z(2)
                            sig_z=dataz_hilb(sigwin_z(1):sigwin_z(2));
                            sig_n=data_n(sigwin_n(1):sigwin_n(2));
                            sig_e=data_e(sigwin_e(1):sigwin_e(2));
                        else
                            continue
                        end
                        % only use events with similar amplitude magnitude
                        if sum(sig_z.^2)/sum(sig_n.^2) >= 0.8 && sum(sig_z.^2)/sum(sig_e.^2) >= 0.8
                            sig_z=sig_z/max(sig_z);
                            sig_n=sig_n/max(sig_n);
                            sig_e=sig_e/max(sig_e);
                        else
                            continue
                        end

                        if length(sig_z)==length(sig_n) && length(sig_z)==length(sig_e)
                            [corr,x,corr2,x2] = calcOrient(sig_z,sig_n,sig_e,baz);
                            iee = iee + 1; % effective event + 1

                            % each column of outdata represents the effective earthquake serial number
                            % Each row of outdata represents the identifier of frequency, identifier of
                            % event, back azimuth, cross correlation, x, second cross correlation and x2.
                            outdata(iee,:) = [ifreq, ie, baz, corr, x, corr2, x2, snr_z];
                        else
                            continue
                        end

                    else
                        disp('Low SNR, skip!');
                    end % end snr
                else
                    disp('No file, skip!')
                end % end exist file

            end % end event loop

        end % end frequency loop


        save(outfile,'outdata','-ascii');


    end % end overwrite


    plotOrient(outfile,mincorr,network,station); % network and station are not necessary

    figure(1)
    set(gcf,'PaperPositionMode','manual');
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperOrientation','portrait');
    set(gcf,'PaperPosition',[.05 1.5 8 5]);
    xlabel('Coherence');
    ylabel('Degrees from N'); % Degrees from the direction of N component
    pdffile1 = sprintf('%s/%s_%s_nor.pdf',figoutpath,network,station);
    print('-dpdf',pdffile1);

    figure(2)
    pdffile2 = sprintf('%s/%s_%s_ros.pdf',figoutpath,network,station);
    print('-dpdf',pdffile2);

end % end station loop