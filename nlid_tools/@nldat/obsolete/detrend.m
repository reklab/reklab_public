function Zout=detrend(ZIn)
% Detrend wrapper for nldat objects

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 


Zout=ZIn;
[nsamp,nchan,nreal]=size(ZIn);
for i=1:nchan,
  for j=1:nreal,
    Zout.dataSet(:,i,j)=detrend(ZIn.dataSet(:,i,j))   
  end
end
Zout.comment =[ ZIn.comment '; detrend'];
return


