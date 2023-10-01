function [] = Z_plot(f, Z, L_initial, plot_title)
%DTFT_PLOT Plot the given function for ex. 4
%   Detailed explanation goes here

figure;
for i = 1 : 4
    L = i + L_initial - 1;
    subplot(2, 2, i)
    plot(f, Z(f, L))
    title(plot_title + " (L = " + num2str(L) + ")")
    xlabel("frequency [Hz]")
    ylabel("magnitude")
    grid on
end

end

