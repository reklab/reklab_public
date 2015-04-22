clear;close all;clc
echo on
% This is demo 4 of the SMI toolbox In this demo we will show the application
% Seperable Least Squares(SLS) techniques, in order to reduse the parameter 
% space during optimization. SLS can be used to optimize linear models, in
% continuous time and discrete time, as well as wiener systems.
%
% A critical step in optimization is finding an accurate initial estimate.
% SLS optimization for linear models uses only the matrices A and C. 
% This is advantageous, since the moesp algorithms estimate B and D 
% only in a second step. For the initial estimate we can skip this step.
%
pause % Press any key to continue.
clc
% This demo shows the use of the functions provided by the toolbox.
% First we start with dslslin. This function does the optimization of 
% linear discrete time systems.
% We use the same system as in smidemo2 to illustrate the use of the 
% optimization functions. This is a three degree of freedom 
% mass-spring-damper system (see manual).
% The output y(k) is perturbed by a process noise w(k), and a
% measurement noise v(k)

den=[1 -0.5 0.7];
denn=[1 0.5 0.5];
num=[1  0.5];
u=prbn(300,0.1)-0.5;
v=randn(300,1);
l=1;m=1;
y=dlsim(num,den,u+0.1*std(u)*randn(300,1));
y=y+dlsim(num,den,v)*0.1*std(y);

% We estimate first a complete system using dordpo/destac. This is our
% initial estimate. 
u=detrend(u);y=detrend(y);
[Sn, R]=dordpo(u,y,12);
[Ai,Ci]=destac(R,6);

pause % Press any key to continue.
clc

%In this case we estimate B and D also, but only to
% verify how good the initial estimate is.
[Be,De]=destbd(u,y,Ai,Ci);
yi=dlsim(Ai,Be,Ci,De,u);
vaf(y,yi)

if max(abs(eig(Ai)))>1
  % make Ai stable, just to fit in the optimization.
  Aold=Ai
  Ai=Ai*(1/max(abs(eig(Ai))))*0.95;
end
model=[1 0 0 0];
options=zeros(14,1);
options(1)=1;
options(14)=2000;

[Ao,Bo,Co,Do,xinit]=dslslin(u,y,Ai,Ci,[],model,options);

% Finally the variance accounted for is calculated. Let's hope it is
% over 99%.

ye=dlsim(Ao,Bo,Co,Do,u);
vaf_end=vaf(y,ye)

echo off
if vaf_end>99
  eval('lalala','')
disp('Hallelujah');
end

% It is also possible to optimize a steady state Kalman filter. This works
% almost exactly identical to the previous example. 
% We can estimate an initial guess with destac, and use this guess as a
% starting value for the optimization in dslslin.


% R is already caluclated in the previous example, so we just use
[Ai,Ci]=destac(R,6);
Ki=destk(Ai,Ci,R)
Ki2=destk(Ai,Be,Ci,De,R)

if max(abs(eig(Ai)))>1
  % make Ai stable, just to fit in the optimization.
  Aold=Ai;
  Ai=Ai/max(abs(eig(Ai)))*0.95;
end

% the third flag in model has to be set to 1 now, to indicate that we want
% to estimate the Kalman filter gain also
model=[1 0 1 0];
[AO,BO,CO,DO,KO,xinit]=dslslin(u,y,Ai,Ci,Ki,model,options);


% Before
yki=dlsim(Ai-Ki*Ci,[Be-Ki*De,Ki],Ci,[De,zeros(l,l)],[u,y]);
vaf_before=vaf(y,yki)

% and after
yke=dlsim(AO-KO*CO,[BO-KO*DO,KO],CO,[DO,zeros(l,l)],[u,y]);
vaf_after=vaf(y,yke)



















