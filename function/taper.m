function taper_data = taper(data,tap_width)
% tap_width: width in percent of frequency vector of cosine taper
% tap_width can be 0, then the taper operation will not be executed.
% Yuechu Wu
% 12131066@mail.sustech.edu.cn
% 2022-10-21

if tap_width == 0
    taper_data = data;
else
    taper_data = data.*tukeywin(length(data),tap_width);
end

return