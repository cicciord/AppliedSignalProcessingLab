close all
clear all
clc

Fs = 44100;
T = 1.25;

progression_root_notes = [587.33 440 493.88 369.99 392.00 293.66 392.00 440];
progression_root_notes_types = [0 0 1 1 0 0 0 0]; 

progression = generate_progression(progression_root_notes, progression_root_notes_types, Fs, T);
progression_obj = audioplayer(progression, Fs);
playblocking(progression_obj)
audiowrite("progression.wav", progression, Fs)


% internal functions

function progression = generate_progression(F0, type, Fs, T)
    progression = zeros(1, T*Fs*length(F0));
    for i = 1 : length(F0)
        chord = generate_chord_adsr(F0(i), Fs, T, type(i));
        progression(1+(i-1)*T*Fs : i*T*Fs) = chord;
    end

end


function chord = generate_chord_adsr(F0, Fs, T, type)
    chord = adsr_envelope(generate_chord(F0, Fs, T, type), Fs, T);
end


function note_env = adsr_envelope(note, Fs, T)
    note_env = note .* interp1([0 0.2 0.4 0.75 1.25], [0 1 0.7 0.7 0], 0 : 1/Fs : T - 1/Fs, "pchip");;
end


function pwm_filt = generate_note(F0, Fs, T)
    Dt = 1/Fs;
    t = 0 : Dt : T - Dt;
    sawt = sawtooth(2*pi*F0*t);
    
    lfo = F0/240;
    lfo_sine = sin(2 * pi * lfo * t);
    
    pwm = sawt < lfo_sine;

    Bt = 15*F0 * 2 * pi / Fs;
    Wp = F0 * 2 * pi / Fs;
    Ws = 16*F0 * 2 * pi / Fs;
    
    Wc = (Wp + Ws) / 2;

    N = ceil((6.1 * pi) / Bt);
    n = 0 : N - 1;

    window = bartlett(N);
    window = window';

    M = ceil((N-1)/2);

    h_id = Wc/pi * sinc(Wc/pi * (n - M));
    
    h = window .* h_id;
    h = h/sum(h);
    
    pwm_filt = filter(h, 1, pwm);
end

function chord = generate_chord(F0, Fs, T, type)
    if type == 1
        F3 = F0 * 2^(3/12);
    else
        F3 = F0 * 2^(4/12);
    end

    F5 = F0 * 2^(7/12);

    chord = generate_note(F0, Fs, T) + generate_note(F3, Fs, T) + generate_note(F5, Fs, T);
    chord = chord / max(abs(chord));
end