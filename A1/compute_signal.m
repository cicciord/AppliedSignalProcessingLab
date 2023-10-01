function [t, x, f, X_psd] = compute_signal(Fs, T, fi, phii, A)
%COMPUTE_SIGNAL Computes the signal x and its psd X
%   compute a signal made up by the sum of 3 sine functions multiplied by a
%   rect function. Also compute its power spectral density
%
%   input parameters:
%   Fs   => sampling frequency
%   T    => rect function duration
%   fi   => array of dimension 1x3 containing the frequencies
%   phii => array of dimension 1x3 containing the phases
%   A    => amplitude of rect function
%
%   return values:
%   t     => time axis
%   x     => time domain resulting signal
%   f     => frequency axis
%   X_psd => psd

% axis definition
Ts = 1 / Fs;
t = 0 : Ts : T - Ts;

% sum of sign signals with random phases
b = sin(2 * pi * fi(1) * t + phii(1)) + sin(2 * pi * fi(2) * t + phii(2)) ...
    + sin(2 * pi * fi(3) * t + phii(3));

% equivalent of the multiplication by a rect function
x = A * b;

% frequency domain
Df = 1 / T;
f = -Fs/2 : Df : Fs/2 - Df;

X = fftshift(fft(x));
X_psd = abs(X).^2;
% X_psd = X_psd / (max(X_psd) / (A^2 * T^2));
X_psd = X_psd / length(X_psd)^2;

end

