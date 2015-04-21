 function y = sum (x, DIM);
 if nargin==1,
    y=sum((x.Data));
 else
    y=sum(x.Data,DIM);
 end
return
% nldat/sum
% rek 23 May  2002

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 
