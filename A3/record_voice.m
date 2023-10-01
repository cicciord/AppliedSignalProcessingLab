close all
clear
clc

filename = "my_voice.wav";

if isfile(filename)
    disp("Audio already recordered, plese delete the file if you intend record it again.")
    [data, Fs] = audioread(filename);
    sound(data, Fs);
else

Fs = 8000;
n_bits = 16;
n_ch = 1;
T = 3;

r = audiorecorder(Fs, n_bits, n_ch);
disp("Recording...")
recordblocking(r, T);
disp("Done!")

audiowrite(filename, getaudiodata(r), Fs);

end
