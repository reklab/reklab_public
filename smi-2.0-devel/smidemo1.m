clear;close all;clc
echo on
% This is demo1 of the SMI toolbox 
% This file gives an walkthrough of the possibilities 
% of the SMI toolbox. We will show the
% identification capabilities of PO-MOESP 

% Example 1:
% First we generate data for a second order system with a 
% binary sequence as input. For a realistic example the 
% output signal is contaminated with 10% process noise, that 
% enters the system with the input, and 10% measurement noise,
% that disturbs the measured output directly.

den=[1 -1.69 0.96];
num=[1  0.5];
u=prbn(300,0.1)-0.5;
y=dlsim(num,den,u+0.1*std(u)*randn(300,1));
y=y+0.1*std(y)*randn(300,1);

echo off
figure
subplot(2,1,1)
plot(u)
axis([0,300,-1.1, 1.1])
title('Input signal')
subplot(2,1,2);
plot(y)
title('Output signal')
echo on
% Press any key to continue.
pause 
clc
% we will now try to identify the system from the measured
% input u, and output y.  We assume the following statespace
% model:
%                  x(k+1) = Ax(k) + Bu(k) + w(k)
%                  y(k)   = Cx(k) + Du(k) + v(k)
%                  x(0)   = x0;
%The SMI toolbox identifies this model in three stages:
% 1) First the order of the system, "n" is estimated. This is the
%    dimension of the A matrix.
%
% 2) Then the A, and C matrix are estimated. Here we can
%    also calculate the Kalman gain if we wish to.
%
% 3) In the final step,  B, D and the initial state x0 is
%    estimated
% Press any key to continue.
pause 
clc;
% step 1:
%    In this step, we supply dordpo with the data, and a 
%    dimension parameter "s". The value of "s" needs to be 
%    bigger than the expeced order of the system, for 
%    instance twice as big. 
%    dordpo, returns a datamatrix R, which is used in the
%    next step, and a singular value vector S, from which 
%    we can try to read the order "n" of the system. 


u=detrend(u);y=detrend(y);
s=10;
[S,R]=dordpo(u,y,s);


% For this demo, we used a very high value for "s", 10. This
% gives a better view of the relation between "s" and the
% order of the system. We clearly see the difference between
% the second and the third singular value, the "gap". In the noise free
% case, only "n" singular values are non-zero. In this case we can clearly
% chose "n" to be 2. Press OK after you selected the desired modelorder.
n=orderselect(S);
echo off
clc
if n~=2
  disp(['You choose a different order, but we will set the order back to' ...
	' 2 now'])
  n=2;
end
echo on
% Step 2:
% In this step, we use the function destac. As input we give
% it the data matrix R, that we got from dordpo, and the
% value of "n", that we found from the singular value plot.
% As output, we obtain A and C.


[A,C]=destac(R,n);
% Press any key to continue.
pause
clc
% Step 3: 
% Once A, and C are calculated, we can use them to
% estimate B and D, and the initial state. This is done
% using the functions destbd.

[B,D,x0]=destbd(u,y,A,C,[1 1 1]);

% Press any key to continue.
pause 
clc
% Now that we have found a model for the data, we can see how
% good this model is. We use the model to make estimate the
% output.

ye=dlsim(A,B,C,D,u,x0);
echo off
figure
ynf=dlsim(num,den,u);
subplot(3,1,1);plot(ynf);title('Real noise-free output');
subplot(3,1,2);plot(ye);title('Estimated output');
subplot(3,1,3);plot(ye-y);title('Estimation error');
echo on
% As a figure of merrit, we take the Variance-Accounted-For(%).
% It is an indication of how close the original signal and
% its estimate resamble each other.
% If they are completely equal, the VAF is 100%. In other cases,
% the VAF is lower, down to minus infinity if the error is
% bigger than the original signal itself.

VAF=vaf(ynf,ye)

% Press any key to continue.
pause 
clc
% This is the end of demo 1. In demo 2  we will show the use of
% multiple data batches in the SMI toolbox.
%
echo off
% Bert Haverkamp, April 1997
% copyright (c) 1997 Bert Haverkamp
% LAST MODIFIED:
%    09/04/97 BRJ created 



