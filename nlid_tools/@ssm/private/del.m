function [d_x] = del(x,nDelay)

% this function adds a delay to the input signal 
% x is the input signal
% nDelay - dealy in samples
% d_x is the delayed singal
% n is the length of delay

n = length(x);
d_x = zeros(size(x));
for i=1:n-nDelay
    d_x(nDelay+i) = x(i);
end
end