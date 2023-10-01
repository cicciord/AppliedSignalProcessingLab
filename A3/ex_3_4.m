close all
clear
clear sound
clc

voice_pitch = 112.6761;

a_I = [0.7 0.1 0.1 0.7];
f_I = [1 2 6 16] * voice_pitch;

a_E = [0.7 0.7 1];
f_E = [1 3 6] * voice_pitch;

a_A = [0.1 0.1 0.1 0.7 0.7 0.1 0.1];
f_A = [1 2 3 4 5 6 8] * voice_pitch;

a_O = [0.7 0.3 0.7 0.1];
f_O = [1 2 3 4] * voice_pitch;

a_U = [1 0.7 0.07];
f_U = [1 2 3] * voice_pitch;

p = [];
T = 2;
Fs = 22050;

[A_tone, t] = fnote(a_A, f_A, p, T, Fs);
[E_tone, ~] = fnote(a_E, f_E, p, T, Fs);
[I_tone, ~] = fnote(a_I, f_I, p, T, Fs);
[O_tone, ~] = fnote(a_O, f_O, p, T, Fs);
[U_tone, ~] = fnote(a_U, f_U, p, T, Fs);

tone_final = [A_tone E_tone I_tone O_tone U_tone];
t_final = [t t+T t+2*T t+3*T t+4*T];

tone_obj = audioplayer(tone_final, Fs);
playblocking(tone_obj, Fs)
audiowrite("aeiou_tone.wav", tone_final, Fs)

a_bell = [1 1 1 1 1 1 1 1];
f_bell = [1 2 2.4 3 4.5 5 5.33 6] * 440;
T_bell = 5;
Fs_bell = 44100;

[tone_bell, t_bell] = fnote(a_bell, f_bell, p, T_bell, Fs_bell);

tone_bell_mod = tone_bell .* exp(-t_bell/0.9);
bell_obj = audioplayer(tone_bell_mod, Fs_bell);
playblocking(bell_obj)
audiowrite("bell_tone.wav", tone_bell_mod, Fs_bell)




