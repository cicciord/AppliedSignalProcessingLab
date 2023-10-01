function [] = DTFT_plot(f, X, L_initial, plot_title)
%DTFT_PLOT Plot the absolute value of the given DTFT for ex. 4

figure;
for i = 1 : 4
    L = i + L_initial - 1;
    subplot(2, 2, i)
    plot(f, abs(X(f, L)))
    title(plot_title + " (L = " + num2str(L) + ")")
    xlabel("frequency [Hz]")
    ylabel("magnitude")
    grid on
end

end

