function stairs (x)
% Overlaiud stairs function for nldat
% options not yet implements
[nsamp,nchan,nreal]=size(x);
yval = x.Data;
if isnan(x.DomainValues),
   xval= ((1:nsamp)-1)*x.DomainIncr + x.DomainStart;
else
   xval=x.DomainValues;
end
stairs (xval,yval);

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 
