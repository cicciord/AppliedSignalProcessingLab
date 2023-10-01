function r = my_xcorr(a, b)
%MY_XCORR Summary of this function goes here
%   Detailed explanation goes here

N = length(b);
M = length(a);
b = [flip(b) zeros(1, M - 1)];
B = zeros(N+M-1, M);
for i = 0 : M-1
    B(:, i+1) = circshift(b, i);
end
r = (B * a')';

end

