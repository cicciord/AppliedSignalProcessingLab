function X = my_DFT(x, N)
%MY_DFT Computes the DFT on the x signal of N samples
%   This function compute the DFT applying the analytical formula.
%
%   input parameters:
%   x => signal
%   N => number of samples
%
%   return values:
%   X => DFT of x

X = zeros(1, N);

% implementation with n = 0 : N-1 require N to be exactly the number of
% samples of x, whit other way allow to compute the DFT also in the case N
% does not coincide with the number of samples. This way it is possible to
% control the spectrum resolution
n = 0 : length(x)-1;
for k = 0 : N-1
    e = exp(-1j * 2 * pi * n * k / N);
    X(k+1) = sum(x .* e);
end

X = circshift(X, floor(N/2));

