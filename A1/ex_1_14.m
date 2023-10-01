close all
clear 
clc

% constants definition
T = 20;
f1 = 200;
f2 = 220;
f3 = 230;
Fc = 2e3;
mu = 0;
sigma_sq = 25;

sigma = sqrt(sigma_sq);

N = T * Fc;
Dt = 1 / Fc;
t = 0 : Dt : T - Dt;

% generate signal
w = mu + sigma * randn(1, N);
x = cos(2 * pi * f1 * t) + cos(2 * pi * f2 * t) + cos(2 * pi * f3 * t) + w;

% estimate PSD
NFFT = 128;
window = hamming(NFFT);
n_overlap = 0;

[Sx, f] = pwelch(x, window, n_overlap, NFFT, Fc, 'centered');
Sx = Sx';

figure
plot(f, 10 * log10(abs(Sx)))
title("Estimated PSD")
xlabel("frequency [Hz]")
ylabel("Magnitude [dB]")
grid on

% find optimal values
NFFT_opt = 1000;
window_opt = hamming(NFFT_opt);
n_overlap_opt = 500;

[Sx, f] = pwelch(x, window_opt, n_overlap_opt, NFFT_opt, Fc, 'centered');
Sx = Sx';

figure
plot(f, 10 * log10(abs(Sx)))
title("Estimated PSD Optimized")
xlabel("frequency [Hz]")
ylabel("Magnitude [dB]")
grid on

% simple periodgram
NFFT_simple = N;
window_simple = rectwin(NFFT_simple);
n_overlap_simple = 0;

[Sx_simple, f_simple] = pwelch(x, window_simple, n_overlap_simple, NFFT_simple, Fc, 'centered');
Sx_simple = Sx_simple';

% bartlet periodgram
M = 25;
NFFT_bartlett = (N - mod(N, M)) / M;
window_bartlett = rectwin(NFFT_bartlett);
n_overlap_bartlett = 0;

[Sx_bartlett, f_bartlett] = pwelch(x, window_bartlett, n_overlap_bartlett, NFFT_bartlett, Fc, 'centered');
Sx_bartlett = Sx_bartlett';

% welch periodgram
M = 25;
NFFT_welch = (2 * N - mod(2 * N, M + 1)) / (M + 1);
window_welch = hamming(NFFT_welch);
n_overlap_welch = floor(NFFT_welch / 2);

[Sx_welch, f_welch] = pwelch(x, window_welch, n_overlap_welch, NFFT_welch, Fc, 'centered');
Sx_welch = Sx_welch';

figure
plot(f_simple, 10 * log10(abs(Sx_simple)))
title("Estimated PSD comparison")
xlabel("frequency [Hz]")
ylabel("Magnitude [dB]")
grid on
hold on
plot(f_bartlett, 10 * log10(abs(Sx_bartlett)), LineWidth=2)
hold on
plot(f_welch, 10 * log10(abs(Sx_welch)))
legend("simple", "bartlett", "welch")

% auto-correlation estimator
[Rx, lags] = xcorr(x, 'biased');
[Rx2] = xcorr(x, 'unbiased');
figure
subplot(2, 1, 1)
plot(lags, Rx)
title("Auto-correlation estimator biased")
xlabel("lags")
ylabel("auto-correlation")
grid on
subplot(2, 1, 2)
plot(lags, Rx2)
title("Auto-correlation estimator unbiased")
xlabel("lags")
ylabel("auto-correlation")
grid on

% filter
N_filter = 165;
M_filter = floor(N_filter / 2);
n = (-M_filter : M_filter)';
h = 0.3 * sinc(0.3 * n) .* hamming(N_filter);
y = filter(h, 1, x);

% periodgram for filtered signal
% simple periodgram
[Sy_simple, ~] = pwelch(y, window_simple, n_overlap_simple, NFFT_simple, Fc, 'centered');
Sy_simple = Sy_simple';

% bartlet periodgram
[Sy_bartlett, ~] = pwelch(y, window_bartlett, n_overlap_bartlett, NFFT_bartlett, Fc, 'centered');
Sy_bartlett = Sy_bartlett';

% welch periodgram
[Sy_welch, ~] = pwelch(y, window_welch, n_overlap_welch, NFFT_welch, Fc, 'centered');
Sy_welch = Sy_welch';

figure
plot(f_simple, 10 * log10(abs(Sy_simple)))
title("Estimated PSD after filtering")
xlabel("frequency [Hz]")
ylabel("Magnitude [dB]")
grid on
hold on
plot(f_bartlett, 10 * log10(abs(Sy_bartlett)), LineWidth=2)
hold on
plot(f_welch, 10 * log10(abs(Sy_welch)))
legend("simple", "bartlett", "welch")

% welch with different windowing
window_welch_rect = rectwin(NFFT_welch);
window_welch_hann = hann(NFFT_welch);

[Sy_welch_rect, ~] = pwelch(y, window_welch_rect, n_overlap_welch, NFFT_welch, Fc, 'centered');
Sy_welch_rect = Sy_welch_rect';

[Sy_welch_hann, ~] = pwelch(y, window_welch_hann, n_overlap_welch, NFFT_welch, Fc, 'centered');
Sy_welch_hann = Sy_welch_hann';

figure
plot(f_welch, 10 * log10(abs(Sy_welch_rect)), LineWidth=2)
title("Estimated PSD")
xlabel("frequency [Hz]")
ylabel("Magnitude [dB]")
grid on
hold on
plot(f_welch, 10 * log10(abs(Sy_welch)), LineWidth=2)
hold on
plot(f_welch, 10 * log10(abs(Sy_welch_hann)))
legend("rectangular", "hamming", "hann")


