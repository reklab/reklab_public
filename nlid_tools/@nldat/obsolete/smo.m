function y = smo(x,n, DIM);
% overlaid smooth function for nldat objects

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 


if nargin==2,
  DIM=1;
end
y=x;
xd=x.Data;
yd=xd;
[nsamp,nchan,nreal]=size(x);
switch DIM
  case 1
    
    for ireal=1:nreal,
      for ichan=1:nchan,
	yd(:,ichan,ireal)=smo(xd(:,ichan,ireal),n);
      end
    end
  case 2
    
    for isamp=1:nsamp,
      for ireal=1:nreal,
	yd(isamp,:,ireal)=smo(xd(isamp,:,ireal),n);
      end
    end
  case 3  
    for isamp=1:nsamp,
      for ichan=1:nchan,
	yd(isamp,ichan,:)=smo(xd(isamp,ichan,:),n);
      end
    end
end

y.Data=yd;
set (y,'Comment','Smoothed');
