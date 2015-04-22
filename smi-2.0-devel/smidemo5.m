clear;close all;clc
echo on
% This is demo 5 of the SMI toolbox
% This example shows how you may use the errors-in-variables (EIV) functions 
% of the SMI toolbox to perform closed-loop identification 
% 
% The model structure for the EIV problem is an innovations model:
% 
%     x(k+1)  = Ax(k) + Bu~(k) + w(k)
%     y~(k)   = Cx(k) + Du~(k) + v(k)
%     u(k)    = u~(k) + f(k)
% Only u(k) and y(k) are available for identification. 
%                ______________       
%    u~(k)       |            |         y(k) 
%    ----------->| LTI system |------------>
%          |     |____________|     /|\ 
%    u(k)  |          /|\            |
%    <------           |             |
%         /|\          |             |
%          |           |             |
%          | f(k)      | w(k)        | v(k)
%
% Here f(k), w(k) and v(k) are zero-mean white noise sequences
% independent of the noise-free input u~(k).  The EIV functions can be 
% applied to both open- and closed-loop data. 
 
% Press any key to continue.
pause
clc 
% The plant is a mechanical system consisting of 2 rotating discs
% and an integrator. Its transfer function (order = 5) is 
Pnum = [0 0.98 12.99 18.59 3.30 -0.02]*1e-3; 
Pden = [1 -4.4 8.09 -7.83 4 -0.86]; 

% A controller has been designed based on H_infinity method
Cnum = [0.61 -2.03 2.76 -1.83 0.49];
Cden = [1 -2.65 3.11 -1.75 0.39];

% The process noise is shaped by the following filter 
Fnum = [0 2.89 11.13 2.74]*0.01;
Fden = [1 -2.7 2.61 -0.9];

% Press any key to continue.
pause
clc 
% We first  build a combined  deterministic and stochastic system 
% with as input respectively the plant input, process noise, input 
% measurement noise and output measurement noise. 
% The two ouputs of the combined system are the noisy plant output 
% and the noisy plant input.
% This system is then put in closed loop with the H_infinity controler.

plant=tf(Pnum,Pden,-1);
controler=tf(Cnum,Cden,-1);
noisefilter=tf(Fnum,Fden,-1);
interconnection=ss([],[],[],[1 0 1; 0 1 0],-1); 
combinedsystem1=parallel(plant,noisefilter,[],[],1,1);
combinedsystem=parallel(combinedsystem1,interconnection,1,1,1,2);
closedloopsystem=feedback(combinedsystem,controler,2,1,-1);

% Press any key to continue. 
pause
clc
% With this closed loop system we simulate 1200 samples.
% An external input, a white noise of standard deviation (std) 1, 
% is used to excite the system. The white noise driving the process
% noise filter $F(z)$ has a std 1/9. The input and output noise are
% both of std 0.01. This results in SNRs of approximately 20dB at 
% the plant output and 5dB at the plant input. 
%
% The figure shows the  simulated input and output.
echo off
N = 1200; 
r = randn(N,1);      % external input 
v= 0.01*randn(N,1);  %output measurement noise
f= 0.01*randn(N,1);  %input measurement noise
w=(1/3)*randn(N,1);  % proces noise
Ucl=[w,r,v,f];
Ycl = lsim(closedloopsystem,Ucl);
u = Ycl(:,2);
y = Ycl(:,1); 
figure 
subplot(2,1,1)
plot(u)
title('plant input u')
subplot(2,1,2)
plot(y)
xlabel('sample index')
title('plant output y')
echo on


% Press any key to continue. 
pause
clc
% Now we have the data needed for identification
% The data available for identification is plant input u(k), 
% plant output y(k) and the external excitation r(k). 
% 
% The first step in identification is to find out about the order of the
% system. The function dordeiv estimates the singular values and
% returns them in the vector S from which the order can be found.


[S,R] = dordeiv(u,y,r,20);

% The order can then be visually inspected with the function orderselect.

n=orderselect(S);
% press the OK button to continu.
echo off
if n~=7
  disp(['A different order was choosen, but we will set the order to 7 now'])
  n=7;
end
echo on


% The singular value plot indicates that the order of the system
% is 7(sometimes it doesn't,sorry but this is a live demo!). 
% The expected order is 8 which is the  sum of the orders of the 
% plant and the process noise shaping filter. 
% One pole cannot be recovered from identification and this pole 
% belongs to the shaping filter. The information contents on the 
% process noise shaping filter is usually very little in 
% closed-loop as the controller is designed to reject this signal.


pause % Press any key to continue. 
clc
% We proceed to estimate the system matrices.
% First the A and C matrices 

[A,C] = destac(R,7);
% For the estimation of B and D we use the reference input as 
% instrumental variable, in order to arive at a consistent estimate.

[B,D] = destbd(u,y,A,C,[],[],r);
% Figure 3 compares the frequency responses of the
% estimated model to the one of the true systems.
echo off
estimatedsystem=ss(A,B,C,D,-1);
[mt,pt,w] = bode(combinedsystem);
[me,pe] = bode(estimatedsystem,w);
figure
loglog(w,squeeze(mt(1,2,:)),w,squeeze(me(1,1,:)));
xlabel('frequency')
legend('True FRF','Estimated FRF')

% Press any key to continue.
pause 
clc
% This is the end of demo 5. We hope you will now be able to use the 
% SMI toolbox.
%
echo off
% Bert Haverkamp, November 1999
% copyright (c) 1997 Bert Haverkamp, Tung Chou





