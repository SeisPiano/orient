function preproc_data = datapreproc(raw_data,dt,varargin)
% Data preprocessing
% usage: dataout = datapreproc(datain,sampling interval,'filtband',[0.02 0.05],'tap_width',0)
% usage: dataout = datapreproc(datain,sampling interval,'filtband',[0.01 0.1],'rmean',false,'rtrend',false)
% perform rmean, rtrend and taper by default
% tap_width can be 0, then the taper operation will not be executed.
%
% Yuechu Wu
% 12131066@mail.sustech.edu.cn
% 2022-09-29
%
% added default parameters
% 2023-03-02, Yuechu Wu


% set up input parser
p = inputParser;
addParameter(p,'filtband',0);
addParameter(p,'tap_width',0.01);
addParameter(p,'rmean',true);
addParameter(p,'rtrend',true);
parse(p,varargin{:});
filtband = p.Results.filtband;
tap_width = p.Results.tap_width;
rmean = p.Results.rmean;
rtrend = p.Results.rtrend;


raw_data(isnan(raw_data)) = 0; % change NaN to 0

if rmean
    rmean_data = raw_data - mean(raw_data);
else
    rmean_data = raw_data;
end

if rtrend
    rtrend_data = detrend(rmean_data);
else
    rtrend_data = rmean_data;
end


if filtband
    filt_data = bpfilt(rtrend_data,dt,filtband(1),filtband(2));
else
    filt_data = rtrend_data;
end

if tap_width
    taper_data = filt_data.*tukeywin(length(filt_data),tap_width);
    preproc_data = taper_data;
else
    preproc_data = filt_data;
end

return