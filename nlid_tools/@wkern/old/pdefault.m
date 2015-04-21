function IOut = pdefault(IIn);
% set default parameters for wkern

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

p=get(IIn,'parameters');
j=size(p);
p(j+1)=param('name','Method','default','Toeplitz','help','Identification method' ,...
   'type','select','limits', {'Toeplitz' 'LS'});
IOut=IIn;
set(IOut,'Parameters',p);
