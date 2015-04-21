function y = emean (x);
% emean - ensemble mean for nldat objects

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

[nsam, nchan, nresp]= size(x);
y=nldat(x);
xm=[];
for i = 1:nchan
   xc=squeeze(x.Data(:,i,:))';
   xm= cat(2,xm, mean(xc)');
end
y.Data=xm;;
set (y,'DomainIncr',1,'Comment','Ensemble Mean');
return
