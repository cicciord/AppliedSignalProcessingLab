close all
clear all
clc

Fs = 44100;
T = 3;

% notes
A_note = 220;
E_note = A_note * 2^(-5/12);
C_note = A_note * 2^(-9/12);
F_note = A_note * 2^(-4/12);
G_note = A_note * 2^(-2/12);
D_note = A_note * 2^(-7/12);
B_note = A_note * 2^(2/12);

% chord
Am = generate_chord_adsr([A_note E_note C_note A_note/4], Fs, T);
F = generate_chord_adsr([A_note F_note C_note F_note/4], Fs, T);
CG = generate_chord_adsr([G_note E_note C_note G_note/4], Fs, T);
Gs = generate_chord_adsr([G_note D_note C_note G_note/4], Fs, T);
G = generate_chord_adsr([G_note D_note B_note/2 G_note/4], Fs, T);
C = generate_chord_adsr([G_note E_note C_note C_note/4], Fs, T);

chords = zeros(1, 28*Fs*T);
chords = add_tone(chords, 0, Am, Fs, T);
chords = add_tone(chords, 1, F, Fs, T);
chords = add_tone(chords, 2, CG, Fs, T);
chords = add_tone(chords, 3, Gs, Fs, T);
chords = add_tone(chords, 3.5, G, Fs, T);
chords = add_tone(chords, 4, Am, Fs, T);
chords = add_tone(chords, 5, F, Fs, T);
chords = add_tone(chords, 6, CG, Fs, T);
chords = add_tone(chords, 7, Gs, Fs, T);
chords = add_tone(chords, 7.5, G, Fs, T);
chords = add_tone(chords, 8, Am, Fs, T);
chords = add_tone(chords, 9, F, Fs, T);
chords = add_tone(chords, 10, C, Fs, T);
chords = add_tone(chords, 11, G, Fs, T);
chords = add_tone(chords, 12, Am, Fs, T);
chords = add_tone(chords, 13, F, Fs, T);
chords = add_tone(chords, 14, C, Fs, T);
chords = add_tone(chords, 15, G, Fs, T);
chords = add_tone(chords, 16, Am, Fs, T);
chords = add_tone(chords, 17, F, Fs, T);
chords = add_tone(chords, 18, C, Fs, T);
chords = add_tone(chords, 19, G, Fs, T);
chords = add_tone(chords, 20, Am, Fs, T);
chords = add_tone(chords, 21, F, Fs, T);
chords = add_tone(chords, 22, C, Fs, T);
chords = add_tone(chords, 23, G, Fs, T);
chords = add_tone(chords, 24, Am, Fs, T);
chords = add_tone(chords, 25, F, Fs, T);
chords = add_tone(chords, 26, C, Fs, T);
chords = add_tone(chords, 27, G, Fs, T);

% melody
A_note = 2*A_note;
E_note = A_note * 2^(-5/12);
C_note = A_note * 2^(-9/12);
F_note = A_note * 2^(-4/12);
G_note = A_note * 2^(-2/12);
D_note = A_note * 2^(-7/12);
B_note = A_note * 2^(2/12);

d = generate_chord_adsr(D_note, Fs, T);
d_2 = generate_chord_adsr(D_note, Fs, T/2);
c = generate_chord_adsr(C_note, Fs, T);
c_2 = generate_chord_adsr(C_note, Fs, T/2);
f_2 = generate_chord_adsr(F_note, Fs, T/2);
f = generate_chord_adsr(F_note, Fs, T);
e = generate_chord_adsr(E_note, Fs, T);
e_2 = generate_chord_adsr(E_note, Fs, T/2);
g_2 = generate_chord_adsr(G_note/2, Fs, T/2);
g = generate_chord_adsr(G_note/2, Fs, T);
gs = generate_chord_adsr(G_note/4, Fs, T);
a = generate_chord_adsr(A_note/2, Fs, T);
as = generate_chord_adsr(A_note/4, Fs, T);
a_2 = generate_chord_adsr(A_note/2, Fs, T/2);
b_2 = generate_chord_adsr(B_note/2, Fs, T/2);

melody = zeros(1, 28*Fs*T);

melody = add_tone(melody, 8 + 0/12, d_2, Fs, T);
melody = add_tone(melody, 8 + 1/12, c, Fs, T);
melody = add_tone(melody, 8 + 3/12, c, Fs, T);
melody = add_tone(melody, 8 + 5/12, c_2, Fs, T);
melody = add_tone(melody, 8 + 6/12, d_2, Fs, T);
melody = add_tone(melody, 8 + 7/12, c, Fs, T);
melody = add_tone(melody, 8 + 9/12, c, Fs, T);
melody = add_tone(melody, 8 + 11/12, c_2, Fs, T);
melody = add_tone(melody, 9 + 0/12, d_2, Fs, T);
melody = add_tone(melody, 9 + 1/12, c, Fs, T);
melody = add_tone(melody, 9 + 3/12, c, Fs, T);
melody = add_tone(melody, 9 + 5/12, c_2, Fs, T);
melody = add_tone(melody, 9 + 6/12, d_2, Fs, T);
melody = add_tone(melody, 9 + 7/12, c, Fs, T);
melody = add_tone(melody, 9 + 9/12, c, Fs, T);
melody = add_tone(melody, 9 + 11/12, c_2, Fs, T);
melody = add_tone(melody, 10 + 0/12, f_2, Fs, T);
melody = add_tone(melody, 10 + 1/12, e, Fs, T);
melody = add_tone(melody, 10 + 3/12, e, Fs, T);
melody = add_tone(melody, 10 + 5/12, e_2, Fs, T);
melody = add_tone(melody, 10 + 6/12, f_2, Fs, T);
melody = add_tone(melody, 10 + 7/12, e, Fs, T);
melody = add_tone(melody, 10 + 9/12, e, Fs, T);
melody = add_tone(melody, 10 + 11/12, g_2, Fs, T);
melody = add_tone(melody, 11 + 0/12, e_2, Fs, T);
melody = add_tone(melody, 11 + 1/12, d, Fs, T);
melody = add_tone(melody, 11 + 3/12, d, Fs, T);
melody = add_tone(melody, 11 + 5/12, g_2, Fs, T);
melody = add_tone(melody, 11 + 6/12, f_2, Fs, T);
melody = add_tone(melody, 11 + 7/12, e, Fs, T);
melody = add_tone(melody, 11 + 9/12, d, Fs, T);
melody = add_tone(melody, 11 + 11/12, g_2, Fs, T);
melody = add_tone(melody, 12 + 0/12, d_2, Fs, T);
melody = add_tone(melody, 12 + 1/12, c, Fs, T);
melody = add_tone(melody, 12 + 3/12, c, Fs, T);
melody = add_tone(melody, 12 + 5/12, c_2, Fs, T);
melody = add_tone(melody, 12 + 6/12, d_2, Fs, T);
melody = add_tone(melody, 12 + 7/12, c, Fs, T);
melody = add_tone(melody, 12 + 9/12, c, Fs, T);
melody = add_tone(melody, 12 + 11/12, c_2, Fs, T);
melody = add_tone(melody, 13 + 0/12, b_2, Fs, T);
melody = add_tone(melody, 13 + 1/12, a, Fs, T);
melody = add_tone(melody, 13 + 3/12, a, Fs, T);
melody = add_tone(melody, 13 + 5/12, a_2, Fs, T);
melody = add_tone(melody, 13 + 6/12, b_2, Fs, T);
melody = add_tone(melody, 13 + 7/12, a, Fs, T);
melody = add_tone(melody, 13 + 9/12, a, Fs, T);
melody = add_tone(melody, 13 + 11/12, a_2, Fs, T);
melody = add_tone(melody, 14 + 0/12, g_2, Fs, T);
melody = add_tone(melody, 14 + 1/12, e, Fs, T);
melody = add_tone(melody, 14 + 3/12, e, Fs, T);
melody = add_tone(melody, 14 + 5/12, e_2, Fs, T);
melody = add_tone(melody, 14 + 6/12, f_2, Fs, T);
melody = add_tone(melody, 14 + 7/12, e, Fs, T);
melody = add_tone(melody, 14 + 9/12, e, Fs, T);
melody = add_tone(melody, 14 + 11/12, g_2, Fs, T);
melody = add_tone(melody, 15 + 0/12, e_2, Fs, T);
melody = add_tone(melody, 15 + 1/12, d, Fs, T);
melody = add_tone(melody, 15 + 3/12, d, Fs, T);
melody = add_tone(melody, 15 + 5/12, g_2, Fs, T);
melody = add_tone(melody, 15 + 6/12, f_2, Fs, T);
melody = add_tone(melody, 15 + 7/12, e, Fs, T);
melody = add_tone(melody, 15 + 9/12, d, Fs, T);
melody = add_tone(melody, 15 + 11/12, g_2, Fs, T);
melody = add_tone(melody, 16 + 0/12, d_2, Fs, T);
melody = add_tone(melody, 16 + 1/12, c, Fs, T);
melody = add_tone(melody, 16 + 3/12, c, Fs, T);
melody = add_tone(melody, 16 + 5/12, c_2, Fs, T);
melody = add_tone(melody, 16 + 6/12, d_2, Fs, T);
melody = add_tone(melody, 16 + 7/12, c, Fs, T);
melody = add_tone(melody, 16 + 9/12, c, Fs, T);
melody = add_tone(melody, 16 + 11/12, c_2, Fs, T);
melody = add_tone(melody, 17 + 0/12, d_2, Fs, T);
melody = add_tone(melody, 17 + 1/12, c, Fs, T);
melody = add_tone(melody, 17 + 3/12, c, Fs, T);
melody = add_tone(melody, 17 + 5/12, c_2, Fs, T);
melody = add_tone(melody, 17 + 6/12, d_2, Fs, T);
melody = add_tone(melody, 17 + 7/12, c, Fs, T);
melody = add_tone(melody, 17 + 9/12, c, Fs, T);
melody = add_tone(melody, 17 + 11/12, c_2, Fs, T);
melody = add_tone(melody, 18 + 0/12, f_2, Fs, T);
melody = add_tone(melody, 18 + 1/12, e, Fs, T);
melody = add_tone(melody, 18 + 3/12, e, Fs, T);
melody = add_tone(melody, 18 + 5/12, e_2, Fs, T);
melody = add_tone(melody, 18 + 6/12, f_2, Fs, T);
melody = add_tone(melody, 18 + 7/12, e, Fs, T);
melody = add_tone(melody, 18 + 9/12, e, Fs, T);
melody = add_tone(melody, 18 + 11/12, g_2, Fs, T);
melody = add_tone(melody, 19 + 0/12, e_2, Fs, T);
melody = add_tone(melody, 19 + 1/12, d, Fs, T);
melody = add_tone(melody, 19 + 3/12, d, Fs, T);
melody = add_tone(melody, 19 + 5/12, g_2, Fs, T);
melody = add_tone(melody, 19 + 6/12, f_2, Fs, T);
melody = add_tone(melody, 19 + 7/12, e, Fs, T);
melody = add_tone(melody, 19 + 9/12, d, Fs, T);
melody = add_tone(melody, 19 + 11/12, g_2, Fs, T);
melody = add_tone(melody, 20 + 0/12, d_2, Fs, T);
melody = add_tone(melody, 20 + 1/12, c, Fs, T);
melody = add_tone(melody, 20 + 3/12, c, Fs, T);
melody = add_tone(melody, 20 + 5/12, c_2, Fs, T);
melody = add_tone(melody, 20 + 6/12, d_2, Fs, T);
melody = add_tone(melody, 20 + 7/12, c, Fs, T);
melody = add_tone(melody, 20 + 9/12, c, Fs, T);
melody = add_tone(melody, 20 + 11/12, c_2, Fs, T);
melody = add_tone(melody, 21 + 0/12, b_2, Fs, T);
melody = add_tone(melody, 21 + 1/12, a, Fs, T);
melody = add_tone(melody, 21 + 3/12, a, Fs, T);
melody = add_tone(melody, 21 + 5/12, a_2, Fs, T);
melody = add_tone(melody, 21 + 6/12, b_2, Fs, T);
melody = add_tone(melody, 21 + 7/12, a, Fs, T);
melody = add_tone(melody, 21 + 9/12, a, Fs, T);
melody = add_tone(melody, 21 + 11/12, a_2, Fs, T);
melody = add_tone(melody, 22 + 0/12, g_2, Fs, T);
melody = add_tone(melody, 22 + 1/12, e, Fs, T);
melody = add_tone(melody, 22 + 3/12, e, Fs, T);
melody = add_tone(melody, 22 + 5/12, e_2, Fs, T);
melody = add_tone(melody, 22 + 6/12, f_2, Fs, T);
melody = add_tone(melody, 22 + 7/12, e, Fs, T);
melody = add_tone(melody, 22 + 9/12, e, Fs, T);
melody = add_tone(melody, 22 + 11/12, g_2, Fs, T);
melody = add_tone(melody, 23 + 0/12, e_2, Fs, T);
melody = add_tone(melody, 23 + 1/12, d, Fs, T);
melody = add_tone(melody, 23 + 3/12, d, Fs, T);
melody = add_tone(melody, 23 + 5/12, g_2, Fs, T);
melody = add_tone(melody, 23 + 6/12, f_2, Fs, T);
melody = add_tone(melody, 23 + 7/12, e, Fs, T);
melody = add_tone(melody, 23 + 9/12, d, Fs, T);

final = chords + 0.6*melody;
final = final/max(abs(final));

smooth_duration = 5;
smooth_end = ones(1, length(final));
smooth_slope = 0 : 1/(Fs*smooth_duration) : 1 - 1/(Fs*smooth_duration);
smooth_slope = [zeros(1, length(smooth_end)-length(smooth_slope)) smooth_slope];
smooth_end = smooth_end - smooth_slope;

final_smooth = final .* smooth_end;
final_obj = audioplayer(final_smooth, Fs);
playblocking(final_obj)
audiowrite("nuvole_bianche.wav", final_smooth, Fs)


% internal functions

function sound = add_tone(sound, start, tone, Fs, T)
    sound = [sound(1:start*Fs*T) sound(start*Fs*T+1:start*Fs*T+length(tone))+tone sound(start*Fs*T+length(tone)+1:end)];
end

function chord = generate_chord_adsr(F0, Fs, T)
    chord = adsr_envelope(generate_chord(F0, Fs, T), Fs, T);
end

function note_env = adsr_envelope(note, Fs, T)

    adsr_envelope_edges = [0 0.01 T-0.09 T];
    adsr_envelope_value = [0 1 0.1 0];
    t = 0 : 1/Fs : T - 1/Fs;
    adsr = (t < adsr_envelope_edges(2)) .* (adsr_envelope_value(2)*t/adsr_envelope_edges(2)) + ...
        (t > adsr_envelope_edges(3)) .* ( -adsr_envelope_value(3) * (t-adsr_envelope_edges(3))/(adsr_envelope_edges(4)-adsr_envelope_edges(3))+adsr_envelope_value(3)) + ...
        (t >= adsr_envelope_edges(2) & t <= adsr_envelope_edges(3)) .* exp((log(0.1))/(T-0.015)*(t-0.01));
    note_env = note .* adsr;
end

function chord = generate_chord(F0, Fs, T)
    chord = zeros(1, T*Fs);
    for i = 1 : length(F0)
        chord = chord + generate_note(F0(i), Fs, T);
    end
    chord = chord / max(abs(chord));
end

function pwm_filt = generate_note(F0, Fs, T)
    Dt = 1/Fs;
    t = 0 : Dt : T - Dt;
    init_tone = (sin(2*pi*F0*t) ...
        + sin(2*pi*F0*t + 0.3/F0) + sin(2*pi*F0*t - 0.3/F0) ...
        + sin(2*pi*F0*t + 0.5/F0) + sin(2*pi*F0*t - 0.5/F0) ...
        + sawtooth(2*pi*F0*t + 0.5/F0) + sawtooth(2*pi*F0*t - 0.5/F0) ...
        + sawtooth(2*pi*F0*t + 0.7/F0) + sawtooth(2*pi*F0*t - 0.7/F0) ...
        );
    
    lfo = F0/16;
    lfo_sine = sin(2 * pi * lfo * t) .* (-t/(20*T));

    lfo_2 = F0/8;
    lfo_sine_2 = 0.12 * sin(2 * pi * lfo_2 * t) + 1;
    
    pwm = init_tone < lfo_sine;
    pwm = pwm .* lfo_sine_2;

    Bt_lp = 7*F0 * 2 * pi / Fs;
    Wp_lp = F0 * 2 * pi / Fs;
    Ws_lp = 8*F0 * 2 * pi / Fs;
    
    Wc_lp = (Wp_lp + Ws_lp) / 2;

    N_lp = ceil((6.1 * pi) / Bt_lp);
    n_lp = 0 : N_lp - 1;

    window_lp = bartlett(N_lp);
    window_lp = window_lp';

    M_lp = ceil((N_lp-1)/2);

    h_id_lp = Wc_lp/pi * sinc(Wc_lp/pi * (n_lp - M_lp));
    
    h_lp = window_lp .* h_id_lp;
    h_lp = h_lp/sum(h_lp);

    Bt_hp = 100 * 2 * pi / Fs;
    Wp_hp = 60 * 2 * pi / Fs;
    Ws_hp = 20 * 2 * pi / Fs;
    
    Wc_hp = (Wp_hp + Ws_hp) / 2;

    N_hp = ceil((6.1 * pi) / Bt_hp);
    n_hp = 0 : N_hp - 1;

    window_hp = bartlett(N_hp);
    window_hp = window_hp';

    M_hp = ceil((N_hp-1)/2);
    
    delta_hp = zeros(1, N_hp);
    delta_hp(M_hp+1) = 1;
    h_id_hp = delta_hp - (Wc_hp/pi * sinc(Wc_hp/pi * (n_hp - M_hp)));
    
    h_hp = window_hp .* h_id_hp;
    h_hp = h_hp/sum(h_hp);
    
    pwm_filt = filter(h_lp, 1, pwm);
    pwm_filt = filter(h_hp, 1, pwm_filt);
end