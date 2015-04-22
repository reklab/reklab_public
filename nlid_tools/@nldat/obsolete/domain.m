function d = domain(x);
% returns domain vlaues for nldat object

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

[nsamp,nchan,nreal]=size(x);
if isnan(x.domainValues),
   d= ((1:nsamp)-1)*x.domainIncr + x.domainStart;
else
   d=x.domainValues;
end
d=d(:);
