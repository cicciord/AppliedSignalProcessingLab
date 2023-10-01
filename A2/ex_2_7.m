close all
clear
clc

Wp = 2;
Ws = 6;
Rp = 1;
As = 70;

Wp_rad = Wp * 2 * pi;
Ws_rad = Ws * 2 * pi;

N = ceil(log10((10^(Rp/10) - 1) / (10^(As/10) - 1)) / (2 * log10(Wp_rad / Ws_rad)));
Wcp_rad = Wp_rad / (10^(Rp/10) - 1)^(1/(2*N));
Wcs_rad = Ws_rad / (10^(As/10) - 1)^(1/(2*N));
Wc_rad = (Wcp_rad + Wcs_rad) / 2;
Wc = round(Wc_rad / (2*pi), 3);

Fs_analog = 40;
w = 0 : 0.1 : Fs_analog - 0.1;
H_analog = 1 ./ (1 + (w/Wc).^(2*N));

figure
plot(w, 20*log10(abs(H_analog)));
title("Square magnitude of analog filter")
xlabel("frequency [Hz]")
ylabel("$\vert H_a(j\Omega)\vert^2_{dB}$", Interpreter="latex")
grid on

p = zeros(1, N/2);
for k = 0 : N/2 - 1
    p(k+1) = Wc_rad * exp(1j * k * pi / N) * exp(1j * pi * ((N+1) / 2) / N);
end

p = [p conj(p)];

[b_analog, a_analog] = zp2tf([], p, Wc_rad^N);

figure
zplane(b_analog, a_analog)
grid on


Fs = 80;

[b_impinvar, a_impinvar] = impinvar(b_analog, a_analog, Fs);
[b_bilinear, a_bilinear] = bilinear(b_analog, a_analog, Fs);

[H_impinvar, f_discrete] = freqz(b_impinvar, a_impinvar, 1024, Fs);
[H_bilinear, ~] = freqz(b_bilinear, a_bilinear, 1024, Fs);
W = 0 : pi*Fs/1024 : pi*Fs - pi*Fs/1024;
[H_analog, w] = freqs(b_analog, a_analog, W);


figure
plot(f_discrete, 20*log10(abs(H_impinvar).^2))
title("Frequency response square magnitude comparison")
xlabel("frequency [Hz]")
ylabel("magnitude^2 [dB]")
grid on
hold on
plot(f_discrete, 20*log10(abs(H_bilinear).^2))
hold on
plot(w/(2*pi), 20*log10(abs(H_analog).^2))
legend("impulse invariance", "bilinear transform", "analog")

try
    [b_hp, a_hp] = iirlp2hp(b_bilinear, a_bilinear, 2*Wc/Fs, 2*36/Fs);
    figure
    freqz(b_hp, a_hp, 1024, Fs);
    
    [b_bp, a_bp] = iirlp2bp(b_bilinear, a_bilinear, 2*Wc/Fs, 2*[16 24]/Fs);
    figure
    freqz(b_bp, a_bp, 1024, Fs);
catch error
    disp("The filter is sampled at a frequency too low to perform the last part of the script")
    disp(error.message)
    disp("Try running the script with Fs equal to 80")
end
