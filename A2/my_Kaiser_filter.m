function [h, N, b] = my_Kaiser_filter(As, B, Fs, Fc, type)
%MY_KAISER_FILTER Summary of this function goes here
%   Detailed explanation goes here

B_rad_norm = B * 2 * pi / Fs;
Fc_norm = Fc / Fs;

if As < 0
    error("As should be >= 0")
end
if As > 50
    b = 0.1102 * (As - 8.7);
elseif As <= 50 && As >= 21
    b = 0.5842 * (As - 21)^0.4 + 0.07886 * (As - 21);
else
    b = 0;
end

if B <= 0
    error("B must be > 0")
end

N = ceil((As - 8) / (2.285 * B_rad_norm));
M = ceil((N-1)/2);
n = 0 : N-1;

w = besseli(0, b .* sqrt(1 - ((2 * n - N + 1) / (N - 1)).^2)) / besseli(0, b);

if strcmp(type, "-lp")
    if length(Fc_norm) ~= 1
        error("Wrong number of cutoff frequencies")
    end

    h_id = 2 * Fc_norm * sinc(2 * Fc_norm * (n-M));
    h = w .* h_id;

elseif strcmp(type, "-hp")
    if length(Fc_norm) ~= 1
        error("Wrong number of cutoff frequencies")
    end
    
    delta = zeros(1, N);
    delta(M+1) = 1;
    h_id = delta - 2 * Fc_norm * sinc(2 * Fc_norm * (n-M));
    h = w .* h_id;

elseif strcmp(type, "-bp")
    if length(Fc_norm) ~= 2
        errordlg("Wrong number of cutoff frequencies")        
    end

    h_id = 2 * Fc_norm(2) * sinc(2 * Fc_norm(2) * (n-M)) - 2 * Fc_norm(1) * sinc(2 * Fc_norm(1) * (n-M));
    h = w .* h_id;

elseif strcmp(type, "-bs")
    if length(Fc_norm) ~= 2
        errordlg("Wrong number of cutoff frequencies")
    end

    delta = zeros(1, N);
    delta(M+1) = 1;
    h_id = delta - 2 * Fc_norm(2) * sinc(2 * Fc_norm(2) * (n-M)) - 2 * Fc_norm(1) * sinc(2 * Fc_norm(1) * (n-M));
    h = w .* h_id;

else
    errordlg("Filter type non existent", "Filter Type")
    return
end

end

