function y=smo(x,numtimes)
% SMO - smooth with a 3-point, zero-phase,  moving average filter 
%
%	y(i) = x(i-1)/4 + x(i)/2 + x(i+1)/4
%
%	USAGE	: y = smo(x,numtimes)
%
%	x	: input vector
%	numtimes: number of times x should be smoothed
%
% if x is a matrix, smo(x,numtimes) will smooth each
% column of x numtimes
% $Revision: 1.5 $

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

% default is a single pass
if nargin < 2
  numtimes = 1;
end

% allocate storage for output
[nr,nc]=size(x);
y = zeros(nr,nc);

% use filtfilt to implement zero phase filter so h  is the "square root" of
% the desired two-sided irf
h = [0.5 0.5]; 

% if multiple passes are needed, do them to the IRF, rather than the data.
if numtimes > 1
  h1 = h;
  for i = 2:numtimes
    h = conv(h,h1);
  end
end
h=h'; 

% we will have to zero pad to fix the ends of the record
% to undo the transient "fixes" provided by filtfilt.
hlen = length(h);
pad = zeros(hlen,1);

for i = 1:nc
  xx = [pad;x(:,i);pad];
  yy = filtfilt(h,1,xx);
  y(:,i) = yy(hlen+1:end-hlen);
end
return
