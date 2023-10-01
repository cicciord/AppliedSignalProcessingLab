function [S, f, t] = my_spectrogram(x, M, L, N, Fs, plot_flag)
%MY_SPECTROGRAM Summary of this function goes here
%   Detailed explanation goes here

if mod(M, 2) == 0
    M = M + 1;
end

window = hamming(M);
window = window';

R = M - L;
Nx = length(x);
T = floor((Nx - L) / (M - L));

f = 0 : Fs/N : Fs/2 - Fs/N;

S = zeros(N/2, T);

t = 0 : R : R*(T-1);

t = t + (M+1)/2;
t = t / Fs;

for i = 0 : T-1
    frame = [zeros(1, (N-M+1)/2) (x(i*R+1 : i*R+M) .* window) zeros(1, (N-M-1)/2)];
    Frame = fft(frame);
    S(:,i+1) = Frame(1:N/2);
end

if plot_flag == 1
    figure
    imagesc(t, f, 20*log10(abs(S)))
    axis xy
end


end

