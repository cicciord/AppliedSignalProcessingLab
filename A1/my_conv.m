function z = my_conv(x, y, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mode = 'l';

if nargin > 3
    error("Too many parameters")
end
if nargin == 3
    param = varargin{1};
    if param ~= 'l' && param ~= 'c'
        error("Wrong parameter, expected 'l' or 'c'")
    end
    mode = param;
end

% linear convolution
% if mode == 'l'     
%     N = length(x);
%     M = length(y);
%     x = [x, zeros(1, N)];
%     y = [y, zeros(1, M)];
%     z = zeros(1, M+N-1);
%     for n = 1 : N+M-1
%         for k = 1 : N
%             if(n-k+1 > 0)
%                 z(n) = z(n) + x(k) * y(n-k+1);
%             end
%         end
%     end

if mode == 'l'
    N = length(x);
    M = length(y);
    x = [x zeros(1, M - 1)];
    B = zeros(N+M-1, M);
    for i = 0 : M-1
        B(:, i+1) = circshift(x, i);
    end
    z = (B * y')';


% circular convolution
elseif mode == 'c'
    N = length(x);
    M = length(y);
    N_z = max(N, M);
    x = [x, zeros(1, N_z - N)];
    y = [y, zeros(1, N_z - M)];

    B = zeros(N_z, N_z);
    for i = 0 : N_z-1
        B(:, i+1) = circshift(y, i);
    end

    z = (B * x')';
end

end

