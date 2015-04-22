 function y = cumsum (x, DIM);

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

y=nldat(x);
if nargin==1,
  y.dataSet=cumsum((x.dataSet));
else
  y.dataSet=cumsum(x.dataSet,DIM);
end

set (y,'comment','Cumsum');
return
% nldat/cumsum
% rek 9 April 2000
