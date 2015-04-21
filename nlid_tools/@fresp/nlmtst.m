function fs=nlmtst(F)
% fresp/nlmtst
%
%fs=nlmtst(fresp) test identification from data
% $Revision: 1.5 $
% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

z=nlid_sim ('L1');
set(z,'comment','Low pass data set');
figure(1);
plot (z); 
fs=fresp(z);
figure(2)
plot(fs);

% Test conversion to nldat
N=nldat(fs);
% test convertsion to IRF
I=irf(fs);
figure(3); 
plot(I);

% Test simulation 

%
r=nlid_resid(fs,z);

%
% Test of delay
fs1=delay (fs,.1);
plot(fs1);
plot(irf(fs1));







% test identification by FFTing IRF;

I = irf(z);
FI=fresp(I);
figure(3); plot (FI);
figure(4); plot(I)



