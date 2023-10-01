close all
clear
clc

filename = 'pulse.txt';
delimiterIn = ' ';
headerlinesIn = 1;
Data_struct = importdata(filename, delimiterIn, headerlinesIn);

Led_R = Data_struct.data(:,1); % RED (R)
Led_IR = Data_struct.data(:,2); % INFRARED (IR)

Fs = 100;
N = length(Led_R);
T = N/Fs;
Dt = 1/Fs;
t = 0 : Dt : 60 - Dt;

r = Led_R(Fs*10:Fs*70-1);
ir = Led_IR(Fs*10:Fs*70-1);

figure
subplot(2, 1, 1)
plot(t, r)
title("RED")
xlabel("time [s]")
subplot(2, 1, 2)
plot(t, ir)
title("INFRARED")
xlabel("time [s]")


% low pass filter
Wp_lp = 3;
Ws_lp = 6;
Rp_lp = 1;
As_lp = 50;

Wp_rad_lp = Wp_lp * 2 * pi;
Ws_rad_lp = Ws_lp * 2 * pi;

N_lp = ceil(log10((10^(Rp_lp/10) - 1) / (10^(As_lp/10) - 1)) / (2 * log10(Wp_rad_lp / Ws_rad_lp)));
Wcp_rad_lp = Wp_rad_lp / (10^(Rp_lp/10) - 1)^(1/(2*N_lp));
Wcs_rad_lp = Ws_rad_lp / (10^(As_lp/10) - 1)^(1/(2*N_lp));
Wc_rad_lp = (Wcp_rad_lp + Wcs_rad_lp) / 2;
Wc_lp = round(Wc_rad_lp / (2*pi), 3);

p = zeros(1, N_lp/2);
for k = 0 : N_lp/2 - 1
    p(k+1) = Wc_rad_lp * exp(1j * k * pi / N_lp) * exp(1j * pi * ((N_lp+1) / 2) / N_lp);
end

p = [p conj(p)];

[b_analog_lp, a_analog_lp] = zp2tf([], p, Wc_rad_lp^N_lp);

Fs = 100;

[b_lp, a_lp] = impinvar(b_analog_lp, a_analog_lp, Fs);

figure
freqz(b_lp, a_lp, 1024, Fs);


r_lp = filtfilt(b_lp, a_lp, r);
ir_lp = filtfilt(b_lp, a_lp, ir);


% high pass
Wp_hp = 0.75;
Wp_rad_hp = Wp_hp * 2 * pi / Fs;
Ws_hp = 0.05;
Ws_rad_hp = Ws_hp * 2 * pi / Fs;

B_hp = (Wp_rad_hp - Ws_rad_hp);

Wc_rad_hp = (Wp_rad_hp + Ws_rad_hp) / 2;
Wc_norm_hp = Wc_rad_hp / (2 * pi);
Wc_hp = Wc_norm_hp * Fs;

N_hp = ceil((6.2 * pi) / B_hp);

n_hp = 0 : N_hp - 1;
window_hp = 0.5 - 0.5 * cos(2 * pi * n_hp / (N_hp-1));

M_hp = ceil((N_hp-1)/2);

delta_hp = zeros(1, N_hp);
delta_hp(M_hp+1) = 1;

h_id_hp = delta_hp - 2 * Wc_norm_hp * sinc(2 * Wc_norm_hp * (n_hp - M_hp));

h_hp = window_hp .* h_id_hp;

figure
freqz(h_hp, 1, 1024, Fs);

r_filt = filtfilt(h_hp, 1, r_lp);
ir_filt = filtfilt(h_hp, 1, ir_lp);

figure
subplot(3, 1, 1)
plot(t, r)
title("original Red signal")
xlabel("time [s]")
subplot(3, 1, 2)
plot(t, r_lp)
title("after low pass")
xlabel("time [s]")
subplot(3, 1, 3)
plot(t, r_filt)
title("after high pass")
xlabel("time [s]")

figure
subplot(3, 1, 1)
plot(t, ir)
title("original Infrared signal")
xlabel("time [s]")
subplot(3, 1, 2)
plot(t, ir_lp)
title("after low pass")
xlabel("time [s]")
subplot(3, 1, 3)
plot(t, ir_filt)
title("after high pass")
xlabel("time [s]")


N_fft = 2^nextpow2(length(r_filt));
R_fft = fftshift(fft(r_filt, N_fft));

Df = Fs/N_fft;
f = -Fs/2 : Df : Fs/2 - Df;

figure
plot(f, 10*log10(abs(R_fft)));
title("Filtered red signal spectrum")
xlabel("frequency [Hz]")
ylabel("magnitude [dB]")
grid on

[~, idx] = max(R_fft);
pulse_rate_hz = abs(f(idx));
pulse_rate_bpm = pulse_rate_hz * 60;
disp("The pulse rate computed is: " + num2str(pulse_rate_bpm) + " BPM")

[r_ac_max, r_ac_max_idx] = findpeaks(r_filt, t);
[r_ac_min, r_ac_min_idx] = findpeaks(-r_filt, t);
r_ac_min = -r_ac_min;
r_ac_max_interp = interp1(r_ac_max_idx, r_ac_max, t, "spline");
r_ac_min_interp = interp1(r_ac_min_idx, r_ac_min, t, "spline");

[r_dc_min, r_dc_min_idx] = findpeaks(-r_lp, t);
r_dc_min = -r_dc_min;
r_dc_min_interp = interp1(r_dc_min_idx, r_dc_min, t, "spline");

[ir_ac_max, ir_ac_max_idx] = findpeaks(ir_filt, t);
[ir_ac_min, ir_ac_min_idx] = findpeaks(-ir_filt, t);
ir_ac_min = -ir_ac_min;
ir_ac_max_interp = interp1(ir_ac_max_idx, ir_ac_max, t, "spline");
ir_ac_min_interp = interp1(ir_ac_min_idx, ir_ac_min, t, "spline");

[ir_dc_min, ir_dc_min_idx] = findpeaks(-ir_lp, t);
ir_dc_min = -ir_dc_min;
ir_dc_min_interp = interp1(ir_dc_min_idx, ir_dc_min, t, "spline");

i_ac_red = r_ac_max_interp - r_ac_min_interp;
i_dc_red = r_dc_min_interp;
i_ac_ired = ir_ac_max_interp - ir_ac_min_interp;
i_dc_ired = ir_dc_min_interp;

R = ( (i_ac_red) ./ (i_dc_red) ) .* ( (i_dc_ired) ./ (i_ac_ired) );
R_avg = mean(R);
SaO2 = 110 - 25 * R_avg;
disp("The saturation level is: " + num2str(SaO2))

figure
subplot(2, 2, 1)
plot(t, r_filt)
hold on
plot(t, r_ac_max_interp, "r")
hold on
plot(r_ac_max_idx, r_ac_max, "ro")
hold on
plot(t, r_ac_min_interp)
hold on
plot(r_ac_min_idx, r_ac_min, "ro")
subplot(2, 2, 2)
plot(t, r_lp)
hold on
plot(t, r_dc_min_interp)
hold on
plot(r_dc_min_idx, r_dc_min, "ro")
subplot(2, 2, 3)
plot(t, ir_filt)
hold on
plot(t, ir_ac_max_interp)
hold on
plot(ir_ac_max_idx, ir_ac_max, "ro")
hold on
plot(t, ir_ac_min_interp)
hold on
plot(ir_ac_min_idx, ir_ac_min, "ro")
subplot(2, 2, 4)
plot(t, ir_lp)
hold on
plot(t, ir_dc_min_interp)
hold on
plot(ir_dc_min_idx, ir_dc_min, "ro")
sgtitle("SaO_2=" + num2str(SaO2) + "%")




