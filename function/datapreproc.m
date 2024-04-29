function data_prep = datapreproc(data,dt,varargin)

% Data preprocessing
% usage: dataout = datapreproc(datain,sampling interval,'filt_band',[0.02 0.05],'tap_width',0)
% usage: dataout = datapreproc(datain,sampling interval,'filt_band',[0.01 0.1],'rmean',false,'rtrend',false)
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
addParameter(p,'filt_band',0);
addParameter(p,'tap_width',0.01);
addParameter(p,'rmean',true);
addParameter(p,'rtrend',true);
parse(p,varargin{:});
filt_band = p.Results.filt_band;
tap_width = p.Results.tap_width;
rmean = p.Results.rmean;
rtrend = p.Results.rtrend;


data(isnan(data)) = 0; % change NaN to 0

if rmean
    data_prep = data - mean(data);
else
    data_prep = data;
end

if rtrend
    data_prep = detrend(data_prep);
end

if tap_width
    data_prep = taper(data_prep,tap_width);
end

if filt_band(1)
    data_prep = bpfilt(data_prep,filt_band(1),filt_band(2),dt);
end


return