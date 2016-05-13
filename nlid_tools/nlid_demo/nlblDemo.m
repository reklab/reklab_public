function i=nlblDemo(i)
% test of NLBL  identification
%% hk - method
clear all
sys = 'N3HP'
disp('hk - test');
[z,m]=nlid_sim (sys);
NHK=nlbl;
set(NHK,'idMethod','hk','nLagLE',32,'maxOrderNLE',3,'displayFlag',true,'threshNSE',.001);
NHK=nlident(NHK,z);
figure(1); plot(NHK); 
figure(2);
nlid_resid(NHK,z);

%% sls method

disp('nlbl identificaton using sls');
NSLS=nlbl;
set(NSLS,'idMethod','sls', 'nIterMax', 10, 'nLagLE', 32,'maxOrderNLE',3,'displayFlag',true, ...
    'nIterMax', 100);
NSLS=nlident(NSLS,z);figure (1)
figure(3);
plot (NSLS) ;
figure(4);
nlid_resid(NSLS,z);;

%% subspace method
disp('nlbl identification using subspace method');
NSS=nlbl;
set(NSS,'idMethod','subspace','nLagLE',32,'maxOrderNLE',3,'displayFlag', true);
NSS=nlident(NSS,z);
figure (5);
plot(NSS)
figure(6);
nlid_resid(NSS,z);

%% Test normalizations
figure(7)
N=normGainLE(NHK);
nlid_resid(N,z); 
plot(N);

figure(8);
N=normCoefLE(NHK);
nlid_resid(N,z);
plot(N);

figure(9);
N=normCoefNLE(NSLS);
nlid_resid(N,z);
plot(N);

%% Test of methods
NS=smo(NHK); plot (NS)
% subref
disp(NHK.idMethod)
% parameter value
disp(NHK.nLagLE)
% test set
disp('nlblDemo done'); 
% nlbl/nlmtst
