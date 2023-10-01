function [S, f] = my_Bartlett(x, N, M)
%MY_BARTLETT Summary of this function goes here
%   Detailed explanation goes here

if M > N
    error('M cannot be larger then N')
end

L = (mod(N, M) ~= 0) * (M - mod(N, M));
N = N + L;

x = [x zeros(1, L)];
D = N / M;

X = zeros(M, D);
for i = 1 : M
    X(i,:) = my_DFT(x(1+(i-1)*D : i*D),D);
end

X = abs(X).^2 / D;

% S = sum(X) / M;
S = mean(X);
S = circshift(S, -1);
f = -1/2 + 1/D : 1/D : 1/2;

end

