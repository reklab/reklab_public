function [nsamp,nchan,nreal]= size (d, DIM)
% overloaded size function for nldatclass.
[nsamp,nchan,nreal]=size(d.dataSet);
if nargout == 0 | nargout==1,
  if nargin ==1,
    nsamp=[ nsamp nchan nreal];
  else
    nsamp=size(d.dataSet, DIM);
  end
elseif nargout ==2
  [nsamp,nchan]=size(d.dataSet);
elseif nargout ==3
  return  
end

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 
