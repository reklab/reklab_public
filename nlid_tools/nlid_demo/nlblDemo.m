function i=nlblDemo(i)
% Demo of NLBL  identification
% Note that to set the paramters of the nonlinear and linear elements you
% must extract the elemnet from the nlbl, change the paramter, and then add
% it back to the nlbl. This is demosntrated below.
%% hk - method
clear all
sys = 'N2L'
disp('hk - test');
[z,m]=nlid_sim (sys);
NHK=nlbl;
set(NHK,'idMethod','hk','displayFlag',true,'threshNSE',.001);
% Set IRF parameters through the second element in the nlbl; 
i=NHK{1,2}; % Get the linear element
set(i,'nLags',50,'irfPseudoInvMode' , 'manual') % Set its parameters; 
NHK{1,2}=i; % Restore the modifed linear element to the nlbl object. 

NHK=nlident(NHK,z);
figure(1); plot(NHK); 
figure(2);
nlid_resid(NHK,z);

%% sls method

disp('nlbl identificaton using sls');
NSLS=nlbl;
set(NSLS,'idMethod','sls', 'nIterMax', 10, 'displayFlag',true, ...
    'nIterMax', 100);

% To set the parameters used to identify the IRF in the IRF elemebnt of
% nlbl.
%
i=NSLS{1,2} % retrieve IRF element
set(i,'irfPseudoInvMode' , 'auto')  % Set the parameters
NSLS{1,2}=i; % Add it back to the nlbl block. 
% SEt the maximum order for the polynomimal 
n=NSLS{1,1}; 
set(n,'polyOrderMax',5)
NSLS{1,1}=n;


NSLS=nlident(NSLS,z);
figure(3);
plot (NSLS) ;
figure(4);
nlid_resid(NSLS,z);;

%% subspace method
disp('nlbl identification using subspace method');
NSS=nlbl;
set(NSS,'idMethod','subspace','displayFlag', true);
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
disp(NHK{1,2}.nLags)
% test set
disp('nlblDemo done'); 
% nlbl/nlmtst
