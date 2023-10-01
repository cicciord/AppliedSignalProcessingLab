function [note, dist] = my_tuner(x, Fs, N)
%MY_TUNER Summary of this function goes here
%   Detailed explanation goes here

A2 = 110;
octave = {'C', -9; 'C^#/D^b', -8; 'D', -7; 'D^#/E^b', -6; 'E', -5; 'F', -4; 'F^#/G^b', -3; 'G', -2; 'G^#/A^b', -1; 'A', 0; 'A^#/B^b', 1; 'B' 2;};

[R, ~] = xcorr(x, 'normalized');

[~, idx] = findpeaks(R(N:end), MinPeakProminence=0.6);
pitch = Fs / (idx(2) - idx(1));

octave_num = 0;
semi_tone = -9;
diff = pitch - A2 * 2^(-2 + -9/12);

for i = -2 : 7
    for j = -9 : 2
        new_diff = pitch - A2 * 2^(i + j/12);
        if abs(new_diff) < abs(diff)
            diff = new_diff;
            octave_num = i+2;
            semi_tone = j;
        else
            break
        end
    end
end

note = [octave{semi_tone + 10}, '_', num2str(octave_num)];
dist = diff;
end

