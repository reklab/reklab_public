function z = min (x,y, DIM);
x=nldat(x);
z=x;
if nargin==1,
  i=find(x.Data==min((x.Data)));
  xmin=x.Data(i);
  d=domain(x);
  DomainValues=d(i);
  set(z,'Data',xmin,'DomainValues',DomainValues);
  
  
elseif nargin==2,
  y=nldat(y);
  z.Data=min(x.Data,y.Data);   
elseif nargin==3,
  y=nldat(y);
  z.Data=min(x.Data,y.Data,DIM);
end

set (z,'Comment','min');
return
% nldat/max

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 
