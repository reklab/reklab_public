function [npath,nel]= size (d)
% overloaded size function for nlm.
% npath - number of parallel paths
% nel - number of elements per path
E=d.Elements;
[npath, nel]=size(E);
if nargout == 0,
   disp(['number of paths =' int2str(npath)]);
   disp(['number of series elements per path =' int2str(nel)]);
   
elseif nargout == 1;
   npath =[ npath nel ];
end
return

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 


