function y = my_filter(b, a, x)
%MY_FILTER Summary of this function goes here
%   Detailed explanation goes here

% check for unstability
if sum(roots(a) >= 1) ~= 0
    errordlg("The filter is unstable", "Filter Error")
    return
end
N = length(x);
Na = length(a);
Nb = length(b);

M = max(Na, Nb) - 1;
x = [zeros(1, M) x];

y = zeros(1, length(x));
af = fliplr(a);
bf = fliplr(b);
for n = 1:N
    y(n+M) = sum(x(n:n+Nb-1) .* bf) - sum(y(n:n+Na-2) .* af(1:end-1));
end
y = y(M+1:end);

