close all
clear
clc

Fs = 400;

% high pass
Fp_hp = 160;
Wp_hp = Fp_hp * 2 * pi / Fs;
Fs_hp = 140;
Ws_hp = Fs_hp * 2 * pi / Fs;
As_hp = 40;

B_hp = (Wp_hp - Ws_hp);

Wc_hp = (Wp_hp + Ws_hp) / 2;
Fc_norm_hp = Wc_hp / (2 * pi);
Fc_hp = Fc_norm_hp * Fs;

N_hp = ceil((6.2 * pi) / B_hp);

n_hp = 0 : N_hp - 1;
window_hp = 0.5 - 0.5 * cos(2 * pi * n_hp / (N_hp-1));

M_hp = ceil((N_hp-1)/2);

delta_hp = zeros(1, N_hp);
delta_hp(M_hp+1) = 1;

h_id_hp = delta_hp - 2 * Fc_norm_hp * sinc(2 * Fc_norm_hp * (n_hp - M_hp));

h_hp = window_hp .* h_id_hp;

[H_hp, f_hp] = freqz(h_hp, 1, 1024, Fs);

figure
plot(f_hp, 20*log10(abs(H_hp)))
title("Frequency response of high-pass filter")
xlabel("frequency [Hz]")
ylabel("magnitude [dB]")
grid on
xline([Fp_hp Fs_hp], ':k', LineWidth=1.6)
xline(Fc_hp, '--r')


% band pass
Fp_bp = [66 94];
Wp_bp = Fp_bp * 2 * pi / Fs;
Fs_bp = [54 106];
Ws_bp = Fs_bp * 2 * pi / Fs;
As_bp = 50;

B_bp = (Wp_bp(1) - Ws_bp(1));

Wc_bp = (Wp_bp + Ws_bp) / 2;
Fc_norm_bp = Wc_bp / (2 * pi);
Fc_bp = Fc_norm_bp * Fs;

N_bp = ceil((6.6 * pi) / B_bp);

n_bp = 0 : N_bp - 1;
window_bp = 0.54 - 0.46 * cos(2 * pi * n_bp / (N_bp-1));

M_bp = ceil((N_bp-1)/2);

h_id_bp = 2 * Fc_norm_bp(2) * sinc(2 * Fc_norm_bp(2) * (n_bp - M_bp)) - 2 * Fc_norm_bp(1) * sinc(2 * Fc_norm_bp(1) * (n_bp - M_bp));

h_bp = window_bp .* h_id_bp;

[H_bp, f_bp] = freqz(h_bp, 1, 1024, Fs);

figure
plot(f_bp, 20*log10(abs(H_bp)))
title("Frequency response of band-pass filter")
xlabel("frequency [Hz]")
ylabel("magnitude [dB]")
grid on
xline([Fp_bp Fs_bp], ':k', LineWidth=1.6)
xline(Fc_bp, '--r')

% stop band
Fp_sb = [102 158];
Wp_sb = Fp_sb * 2 * pi / Fs;
Fs_sb = [118 142];
Ws_sb = Fs_sb * 2 * pi / Fs;
As_sb = 70;

B_sb = (Ws_sb(1) - Wp_sb(1));

Wc_sb = (Wp_sb + Ws_sb) / 2;
Fc_norm_sb = Wc_sb / (2 * pi);
Fc_sb = Fc_norm_sb * Fs;

N_sb = ceil((11 * pi) / B_sb);

n_sb = 0 : N_sb - 1;
window_sb = 0.42 - 0.5 * cos(2 * pi * n_sb / (N_sb-1)) + 0.08 .* cos(4 * pi * n_sb / (N_sb-1));

M_sb = ceil((N_sb-1)/2);

delta_sb = zeros(1, N_sb);
delta_sb(M_sb+1) = 1;

h_id_sb = delta_sb - ((2 * Fc_norm_sb(2) * sinc(2 * Fc_norm_sb(2) * (n_sb - M_sb)) - 2 * Fc_norm_sb(1) * sinc(2 * Fc_norm_sb(1) * (n_sb - M_sb))));

h_sb = window_sb .* h_id_sb;

[H_sb, f_sb] = freqz(h_sb, 1, 1024, Fs);

figure
plot(f_sb, 20*log10(abs(H_sb)))
title("Frequency response of band-stop filter")
xlabel("frequency [Hz]")
ylabel("magnitude [dB]")
grid on
xline([Fp_sb Fs_sb], ':k', LineWidth=1.6)
xline(Fc_sb, '--r')


