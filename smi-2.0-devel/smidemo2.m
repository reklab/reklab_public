 clear;close all;clc
echo on
% Welcome to demo 2 of the SMI toolbox.  In this demo we 
% will present an example of the use  of the SMI toolbox 
% with multiple data sets.

% First we define our system. This is a three degree of
% freedom mass-spring-damper system (see manual).
% The output y(k) is perturbed by a process noise w(k), and a
% measurement noise v(k)

A = [[.9856 .1628;-.1628 .9856] zeros(2,4);
    zeros(2,2) [.8976 .4305;-.4305 .8976] zeros(2,2);
    zeros(2,4) [.8127 .5690;-.5690 .8127]];
B = [.0011;.0134;-.0016;-.0072;.0011;.0034];
C = [1.5119 0 2 0 1.5119 0;
     1.3093 0 0 0 -1.3093 0];
D = zeros(2,1);
Q = 10^(-4)*diag([.0242 3.5920 .0534 1.034 .0226 .2279]);
R = 10^(-2)*diag([2.785 2.785]);
pause % Press any key to continue.
clc
% We have the availability over three experimental data sets for
% identification. The first with a pink input noise u1, the second
% with a multi-sine input u2, and the third with a binary sequence u3.
% A fourth is used for validation

N = 800;T=(0:N-1)';
[b,a]=butter(4,.8);
u1 = filter(b,a,randn(N,1));
w1 = randn(N,6)*sqrt(Q);v1 = randn(N,2)*sqrt(R);
ynf1 = dlsim(A,B,C,D,u1);y1=ynf1+dlsim(A,eye(6),C,zeros(2,6),w1)+v1;
omega2=(0.2:0.4:1.8);
phase2=[.2019,-1.0538,2.5054,-1.101,-.2132];
u2=zeros(N,1);for j=1:5,u2=u2+sin(omega2(j)*T+phase2(j));end
w2 = randn(N,6)*sqrt(Q);v2 = randn(N,2)*sqrt(R);
ynf2 =  dlsim(A,B,C,D,u2);y2=ynf2+dlsim(A,eye(6),C,zeros(2,6),w2)+v2;
u3=-1+2*prbn(N,0.2);w3 = randn(N,6)*sqrt(Q);v3 = randn(N,2)*sqrt(R);
ynf3 =  dlsim(A,B,C,D,u3);y3=ynf3+dlsim(A,eye(6),C,zeros(2,6),w3)+v3;
u4 = filter(b,a,randn(N,1));
w4 = randn(N,6)*sqrt(Q);v4 = randn(N,2)*sqrt(R);
ynf4 = dlsim(A,B,C,D,u1);y4=ynf4+dlsim(A,eye(6),C,zeros(2,6),w1)+v4;
echo off
figure
subplot(3,1,1);plot(u1),title('u1');
subplot(3,1,2);plot(u2);title('u2');
subplot(3,1,3);plot(u3),title('u3');
axis([0,800,-1.1,1.1])
echo on
pause % Press any key to continue.
clc

% We will now use the three data batches to estimate the
% system in three steps:

% First we will build the data matrix R from the three data
% batches, and estimate the order of the system. 
u1=detrend(u1);y1=detrend(y1);
u2=detrend(u2);y2=detrend(y2);
u3=detrend(u3);y3=detrend(y3);
[Sn1, R1]=dordpo(u1,y1,12);
[Sn2, R2]=dordpo(u2,y2,12,R1);
[Sn3, R3]=dordpo(u3,y3,12,R2);
% We plot the singular values, obtained after each step.
% We see that with every added data set, the order of the system
% becomes more clear: 6
echo off
figure
subplot(1,3,1);semilogy(Sn1,'rx');title('Sing. Val. after batch1');
subplot(1,3,2);semilogy(Sn2,'rx');title('Sing. Val. after batch2');
subplot(1,3,3);semilogy(Sn3,'rx');title('Sing. Val. after batch3');
echo on  
pause % Press any key to continue.
clc

% Next we estimate A and C from the last obtained data matrix
% R3
[Ae,Ce]=destac(R3,6);

% With this estimate of A and C, we can now go ahead, and
% try to find B and D. 
% For this, we use the function destbd
[Be,De,R1bd]=destbdx(u1,y1,Ae,Ce,[1 1 0]);
[Be,De,R2bd]=destbdx(u2,y2,Ae,Ce,[1 1 0],R1bd);
[Be,De,x0e]=destbdx(u3,y3,Ae,Ce,[1 1 1],R2bd);

pause % Press any key to continue.
clc

% Now we would like to verify the accuraccy of our model.

% First we plot the poles of the system and the poles of the model in
% One figure.
echo off
figure
angle=(0:pi/50:2*pi);
plot(sin(angle),cos(angle),'k:');
hold on
plot(eig(A),'b+')
plot(eig(Ae),'rx');
axis([-1.1,1.1,-1.1,1.1])
axis('square')

%legend('','Real poles of system','Estimated poles');
echo on
% Because the system is marginally stable, there is a fair chance that
% the identified model is in unstable. Simulating with such an unstable
% model is not really possible.

pause % Press any key to continue.
clc
% Because it seems not possible to predict the ouput on the basis of the
% input alone, we will now estimate a Kalman-filter, to predict the
% output of the system on the basis of the input, and the past output.

Ke=destk(Ae,Ce,R3);

% We generate a new input output set to test our model on:

% Also we calculate the initial state for the validation data batch:
x04e=destx(u4,y4,Ae,Be,Ce,De);
% With these estimates, we calculate the one step ahead prediction of
% the output;
yek = dlsim(Ae-Ke*Ce,[Be-Ke*De Ke],Ce,[De zeros(2,2)],[u4 y4],x04e);
% To evaluate  our estimate, we compare it with the true onestep ahead
% prediction, based on the real system matrices A, B, C, D, Q, R and S.

L = dlqe(A,eye(6),C,Q,R);
K = A*L;
yk = dlsim(A-K*C,[B-K*D K],C,[D zeros(2,2)],[u4 y4]);
vaf(y4,yek) 
vaf(y4,yk)
% This is the end of demo 2. In demo 3 we will show the use of
% the SMI toolbox with a nonlinear system.
%
echo off
% Bert Haverkamp, April 1997
% copyright (c) 1997 Bert Haverkamp
