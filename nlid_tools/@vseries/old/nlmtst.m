function vs=nlmtst(i)
%
% test of vseries identification
%

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

z=nlid_sim ('LN2');
% vseries - fast orthogoal
vs=vseries('nLags',20);
figure(1);
plot(vs);
figure(2)
nlid_resid ( vs, z);

% vkern/nlmtst
