function IOut = pdefault(IIn);
% set default parameters for wkern
p=get(IIn,'parameters');
j=size(p);
IOut=IIn;
set(IOut,'Parameters',p);

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 
