function [wsls , wstp] =nlmtst(i)
%
% test of wseries identification
%

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

z=nlid_sim ('N2L');
% wkernel = Lee-Schetzen
wsls=wseries(z,'nlags',32,'method','LS');
wstp=wseries(z,'nlags',32,'method','Toeplitz');

figure(1);
plot(wsls);
figure(2);
nlid_resid(wsls,z);
figure(3)
plot(wstp);
figure(4)
nlid_resid(wstp,z);

% wseries/nlmtst
