function vk=nlmtst(i)
%
% test of vkern objects%
% vkernal - fast orthogoal
z=rand(10,10);
vk=vkern(z);
plot(vk);
% vkern/nlmtst

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 
