function [tone, t] = fnote(a, f, p, T, Fs)
%FTONE Summary of this function goes here
%   Detailed explanation goes here
if length(a) ~= length(f)
    errordlg("amplitudes and frequencies must have same dimensions", "Error")
    return
end

N = length(a);
M = length(p);
if M < N
    p = [p zeros(1, N-M)];
elseif M > N
    p = p(1:N);
end

Dt = 1/Fs;
t = 0 : Dt : T - Dt;

tone = zeros(1, length(t));
for i = 1 : N
    harmonic = a(i) * cos(2 * pi * f(i) * t + p(i));
    tone = tone + harmonic;
end

tone = tone / max(abs(tone));
end

