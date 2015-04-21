function phi = phixy (x,y,hlen);
%  first-order cross-correlation between two signals 
%
%  Syntax: phi = phixy (x,y,hlen); 
%          phi = phixy (x,hlen);              (autocorrelation)
%
%  Computes the first-order cross-covariance between x and y out to a 
%  memory length of hlen samples.  Computations are performed in the 
%  frequency domain, as in covf -- however unlike that routine, only
%  the cross-correlation is computed.
%

% Copyright 1991-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 

data_len = length(x);

%  if less than three signals are passed, compute the auto-correlation
if nargin < 3
  hlen = y;
  y = x;
 end

x = x - mean(x);
y = y - mean(y);

%  least power of 2 greater than data_len
fft_size = 2.^ceil(log(2*data_len)/log(2));
X=fft(([x(:) y(:); zeros(data_len,2)]),fft_size);

%  calculate the cross-covariance in the frequency domain
XY = X(:,1).*conj(X(:,2));

%  and fft, back to the time domain
xy=fft(XY,fft_size);

%  and take the real part
phi=real(xy(1:hlen))/(data_len*fft_size);

