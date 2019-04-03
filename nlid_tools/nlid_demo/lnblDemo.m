function i=nlblDemo(i)
% test of NLBL  identification
%% Generate test data set
clear all
sys = 'LN3'
disp('LN3 - test');
R=randvar('randvarType','Uniform','minVal',-5,'maxVal',5);
U=nlsim(R,[0:1000]');set(U,'domainIncr',.01);
[z,m]=nlid_sim (sys,U);
N=lnbl;
% set legnth of IRF
I=N{1,1}; set(I,'nLags',50);
N{1,1}=I;

%% Busgang's method
set(N,'idMethod','busgang');

N=nlident(N,z);
figure(1); plot(N); 
figure(2);
nlid_resid(N,z);

%% hk method

set(N,'idMethod','hk');

N=nlident(N,z);
figure(1); plot(N); 
figure(2);
nlid_resid(N,z);


%% phk method

set(N,'idMethod','phk');

N=nlident(N,z);
figure(1); plot(N); 
figure(2);
nlid_resid(N,z);

%% lm method

set(N,'idMethod','lm');

N=nlident(N,z);
figure(1); plot(N); 
figure(2);
nlid_resid(N,z);



% test set
disp('nlblDemo done'); 
% nlbl/nlmtst
