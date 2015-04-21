function out=hard_limit (x,h_min,h_max)
%   Hard limiter static nonllinearity
%
%   Syntax: y = hard_limit (x, h_min, h_max);
%
%   Applies a hard limiter to the input.
%   Bounds of the hard limiter are passed as arguments.
%
%       y   	: hard limiter output
%       x    	: input signal
%       h_min	: minimum value passed by the limiter
%       h_max  	: maximum value passed
%
%   If only one value is passed, the hard limiter is assumed to by
%   symmetric about 0
%
%   No local functions are called
%

% Copyright 1991-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 


if nargin == 2
    h_max = - h_min;
    end
if h_min > h_max
    temp = h_min;
    h_min = h_max;
    h_max = temp;
    end
%
%  cut bottom end
%
x = ((x - h_min) + abs (x - h_min))/ 2 + h_min;
%
%  cut top end
%
out = ((x - h_max) - abs (x - h_max))/ 2 + h_max;


