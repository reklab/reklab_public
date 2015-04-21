function i=nlmtst(i)
%
% test of lnlbl identification
%


% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

clear all
z=nlid_sim ('LNL');
LNL=lnlbl;

LNL=nlident(LNL,z);
figure(1);
plot(LNL); 
nlid_resid(LNL,z);

figure(2)
plot(i)

%error ('nlmtst not yet implemented for lnlbl');

% lnbl/nlmtst
