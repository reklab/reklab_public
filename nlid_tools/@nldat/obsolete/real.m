function z = real(x);
% real  function for nldat variables;
% $Revision: 1.1 $ 

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

z=x;
z.Data=real(x.Data);
