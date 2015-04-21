function cost = mdl_cost(y,yest,d,hlen);
%  MD cost function evaluation
%
%  syntax:  cost = mdl_cost(y,yest,d,hlen);
%
%  d is the number of model parameters

% Copyright 1994-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 


if nargin < 4
  hlen = 0;
 end

y = double(y);
yest = double(yest);
 
N = length(y);
y = y(hlen+1:N);
yest = yest(hlen+1:N);

cost = sum((y-yest).^2)*(1+d*log(N)/N)/N;
