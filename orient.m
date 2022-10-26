% Determination of ocean bottom seismograph (OBS) orientation via 
% Rayleigh-wave polarization.
% 
% original version was written by Dr. Ba Manh Le.
% The input earthquake events must include the information of epicenter
% distance and azimuth in SAC header variables.
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

INPUTdir = 'DATA'; % directory for data input
OUTPUTdir = 'DATA/ORIENT'; % directory for data output
FIGUREdir = 'FIGURES'; % directory for figure output

% If you don't know which horizontal channel corresponds to the N component, 
% try it. Only the correct channel can get reasonable results.
chz = 'HHZ'; %  Z  component channel name
chn = 'HH2'; % N-S component channel name
che = 'HH1'; % E-W component channel name

bpass = [0.015, 0.025, 0.035, 0.045];  % frequency used to calculate
dbp = 0.005;  % the band for each frequency bpass +/- dbp
minsnr = 10;  % signal to noise ratio threshold
velocity = [3, 4.5]; % lower and upper limits of Rayleigh wave velocity
mincorr = 0.8;  % minimum cross correlation use to pick the result
tap_width = 0.1; % width in percent of frequency vector of cosine taper
snr_op = 2; % options for using SNR calculations: (1) sqrt(∑signal^2/∑noise^2) 
                                                % (2) max(abs(signal)/mean(abs(noise))

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

filesuff = sprintf('*%s*%s',chz,ext); % change this for matching with file name format

for ista = 1:length(stations) % begin station loop
    station = stations{ista};
    close all; clear outdata ide;
    filenames = dir(sprintf('%s/%s/%s/%s',INPUTdir,network,station,filesuff));
    outdata=[];
    ide = 0; % effective event

    for ifreq = 1:length(bpass) % begin frequency loop
        lowfreq = bpass(ifreq)-dbp; higfreq = bpass(ifreq)+dbp;

        for ie = 1:length(filenames)  % begin event loop
            filenamez = sprintf('%s/%s/%s',inpath,station,filenames(ie).name);
            filenamen = strrep(filenamez,chz,chn);
            filenamee = strrep(filenamez,chz,che);

            if exist(filenamen,'file') && exist(filenamee,'file')

                fprintf('loading %s\n',filenamez);
                rawdataz = rsac(filenamez);
                [otimez,btimez,dist,baz,az,dt] = lh(rawdataz,'O','B','DIST','BAZ','AZ','DELTA');
                dataz = datapreproc(rawdataz(:,2),dt,[lowfreq,higfreq],tap_width);

                if snr_op == 1
                    snrz = calsnr(dataz,otimez,btimez,dist,dt);
                elseif snr_op == 2
                    snrz = calsnr2(dataz,otimez,btimez,dist,dt);
                else 
                    error('Incorrect snr_op!');
                end
                
                dataz_hilb = imag(hilbert(dataz));

                fprintf('loading %s\n',filenamen);
                rawdatan = rsac(filenamen);
                [otimen,btimen,~,~,~,dt] = lh(rawdatan,'O','B','DIST','BAZ','AZ','DELTA');
                datan = datapreproc(rawdatan(:,2),dt,[lowfreq,higfreq],tap_width);
%                 snrn = calsnr(datan,otimen,btimen,dist,dt);

                fprintf('loading %s\n',filenamee);
                rawdatae = rsac(filenamee);
                [otimee,btimee,~,~,~,dt] = lh(rawdatae,'O','B','DIST','BAZ','AZ','DELTA');
                datae = datapreproc(rawdatae(:,2),dt,[lowfreq,higfreq],tap_width);
%                 snre = calsnr(datae,otimee,btimee,dist,dt);


                if snrz > minsnr
                    sigwinz = round((dist./fliplr(velocity)+otimez-btimez)/dt)+1;
                    sigwinn = round((dist./fliplr(velocity)+otimen-btimen)/dt)+1;
                    sigwine = round((dist./fliplr(velocity)+otimee-btimee)/dt)+1;

                    if length(datae)>=sigwine(2) && length(datan)>=sigwinn(2) && length(dataz_hilb)>=sigwinz(2)
                        sigz=dataz_hilb(sigwinz(1):sigwinz(2));
                        sign=datan(sigwinn(1):sigwinn(2));
                        sige=datae(sigwine(1):sigwine(2));
                    else
                        continue
                    end
                    % only use events with similar amplitude magnitude
                    if sum(sigz.^2)/sum(sign.^2) >= 0.9 && sum(sigz.^2)/sum(sige.^2) >= 0.9
                        sigz=sigz/max(sigz);
                        sign=sign/max(sign);
                        sige=sige/max(sige);
                    else
                        continue
                    end

                    if length(sigz)==length(sign) && length(sigz)==length(sige)
                        [corr,x,corr2,x2] = calcOrient(sigz,sign,sige,baz);
                        ide = ide + 1; % effective event + 1

% each column of outdata represents the effective earthquake serial number
% Each row of outdata represents the identifier of frequency, identifier of 
% event, back azimuth, cross correlation, x, second cross correlation and x2.
                        outdata(ide,:) = [ifreq, ie, baz, corr, x, corr2, x2, snrz];
                    else
                        continue
                    end

                else
                    disp('Low SNR, skip!');
                end
            else
                disp('No file, skip!')
            end

        end % end event loop

    end % end frequency loop

    outfile = sprintf('%s/%s.%s.txt',outpath,network,station);
    save(outfile,'outdata','-ascii');
    pdffile1 = sprintf('%s/%s.%s.nor.pdf',figoutpath,network,station);
    pdffile2 = sprintf('%s/%s.%s.ros.pdf',figoutpath,network,station);
    plotOrient(mincorr,outfile,pdffile1,pdffile2);

end % end station loop