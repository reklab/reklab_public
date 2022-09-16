function [d_x] = del(x,f,s)

% This function adds a delay to an input signal (x) that is not a nldat object
% For inputs of nldat object type, there is a dedicated del method in the nldat class itself. 
% x is the input signal
% d_x is the delayed singal
% f is the sampling frequency
% s is the delay in seconds

%% Comments/Updates:
%++ Ehsan Sobhani: On Feb. 24, 2014, added the support for negative delay values (i.e., advancing a signal in time)
%++ Ehsan Sobhani: On May  09, 2014, I changed the "fix" function for calculating s to "round" function because it was giving me errors sometimes. 

%% Body of the code
n = length(x);
d_x = zeros(size(x));

% s = floor(s*f);
s = round(s*f);

if s>0
    for i=1:n-s
        d_x(s+i) = x(i);
    end
else
    s = -s;
    for i=1:n-s
        d_x(i) = x(i+s);
    end
end
