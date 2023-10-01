close all
clear
clc

% input dialog
prompt = {'Enter A (from 0.1 to 10):', 'Enter T (integer from 1 to 5):'};
dlgtitle = "Rectangular Pulse Values";
dims = [1 50];
definput = {'1', '1'};
answer = inputdlg(prompt, dlgtitle, dims, definput);

A = str2double(answer{1});
T = str2double(answer{2});

% input checks
if A < 0.1 || A > 10
    error('A value is not in the allowed boundaries');
end

if T < 1 || T > 5
    error('T value is not in the allowed boundaries');
end

if not(floor(T) == T)
    error('T is not an integer number');
end

% variables definition
T_max = 10 * T;
Fs = 1000 / T;
Ts = 1 / Fs;

% number of samples
% sampling frequency s.t. there are 1000 samples in the pulse duration
% total duration is 10 times the duration of a pulse
N = Fs * T_max;

% axis definition
t = 0 : Ts : T_max - Ts;

Df = Fs / N;
f = -Fs/2 : Df : Fs/2 - Df;

% signal definition
x = A * rectangularPulse(0, T, t);

% energy of the signal
E_trapz = trapz(t, abs(x).^2);
E_th = A^2 * T;
E_err = abs(E_th - E_trapz) / E_th;

disp(['The energy of the signal computed with the trapezoid method is :    ', num2str(E_trapz)])
disp(['The energy of the signal computed with theoretically is        :    ', num2str(E_th)])
disp(['Relative error in the numerical integral computation           :    ', num2str(E_err)])

% signal plot
figure
plot(t, x, LineWidth=2)
title('Rectangular Pulse')
xlabel('time [s]')
ylabel('x(t)')
grid on
ylim([-0.1*A A*1.1])

% fourier transform of the signal
X = fftshift(fft(x));
X_psd = abs(X).^2;

% amplitude normalization to (A*T)^2
X_psd = X_psd / (max(X_psd) / (A^2 * T^2)); % check if this is right way to do it

% frequency plot
figure
plot(f, X_psd, LineWidth=2);
title('Signal PSD')
xlabel('frequency [Hz]')
ylabel('X(f)')
xlim([-5/T 5/T])
grid on

% energy in frequency domain
Ef = trapz(f, X_psd);
disp(' ')
disp(['The energy of the signal in frequency domain is                :    ', num2str(Ef)]);

Ef_err = abs(E_th - Ef) / E_th;
disp(['Relative error in the numerical integral computation           :    ', num2str(Ef_err)])

% % compute energy at each lobe
n = 1/T : 1/T : 100/T; % lobes number axis
f0 = length(f) / 2 + 1;
step = (1 / T) / Df;
E_n = zeros(1,100);
for i = 1:100
    E_n(i) = trapz(f(f0 + step * (i - 1) : f0 + step * i), X_psd(f0 + step * (i - 1) : f0 + step * i));
end
E_n = cumsum(E_n);

% comparison with total energy (percentage)
E_p = 2 * E_n / Ef; % factor of 2 since only positive side is considered

% lobes energy plots
figure
plot(n(1:10), E_p(1:10), '-o')
title('Energy contained in the first 10 lobes')
xlabel('lobe number')
ylabel('lobe energy over total energy')
grid on

figure
plot(n(10:end), E_p(10:end))
title('Energy contained in lobes from 10 to 100')
xlabel('lobe number')
ylabel('lobe energy over total energy')
grid on

