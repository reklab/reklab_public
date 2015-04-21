function y = filter_ts(fil,x,numsides, incr )
%  performs one and two-sided filtering.
%
%  Syntax: y = filter_ts(fil,x,numsides,incr)
%
% fil       : filter impulse response
% x         : signal to filter
% numsides  : = 1 for a causal filter
%           : = 2 for an anticausal filter (default)
%           ; 
% incr      : samping increment for filter and data (default is 1)
%
%
% Note that for a one-sided filter, the 
% first length(fil) points are garbage, and for a two
% sided filter, the first and last length(fil)/2 
% points are useless.
%
% for one sided filtering, this calls y=filter(fil,1,x);
%

% Copyright 1991, Robert E Kearney, David T Westwick and Eric J Perreault
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

[ri,ci]= size(x);
if (ci > 1)
     	x = x';
end
if nargin < 4
  incr = 1;
end
if nargin < 3
  numsides = 2;
end
numpts = length(x);
halflen = ceil(length(fil)/2);
if numsides == 2
	x = [x ; zeros(halflen,1)];
	y = filter(fil,1,x);
	y = y(halflen:numpts + halflen - 1);
else 
	y=filter(fil,1,x);
end
if (ci > 1)
     	y = y';
end
y = y*incr;
return

