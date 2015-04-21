function v = vaf(x,y)
% VAF - compute variance accounted for between two signals
% v = vaf(x,y)
% using cov.
%
% $Revision: 1.3 $

% Copyright 1994-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

resid=x-y;
v=100*(1-(std(resid)/std(x))^2);
