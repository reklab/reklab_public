clear;close all;clc
echo on
% This is demo3 of the SMI toolbox 
% In this example we show the capabilities of the SMI toolbox
% to identify nonlinear systems, if they are of the Wiener type. 
%  The example is taken from
% T. Wigren 'Convergence analysis of recursive identification algorithms
% based on the nonlinear wiener model'
% [IEEE trans. autom. control. Vol 39 No. 11, 1995]


% A wiener system is  a concatenation of a linear dynamic
% system, and a static nonlinearity.
%       ______________        _______________________
%  u(k) |            |  y(k)  |                     |  z(k)
% ----->| LTI system |------->| static nonlinearity |------>
%       |____________|        |_____________________|
%	
% First we demonstrate the use of the PI-moesp routine to estimate 
% the LTI part of such a system if the input is gausian. The PI-moesp
% algorithm is able to estimate the LTI part unbiased in such a case.

pause % Press any key to continue. 
clc
% system transfer function
num = [.1578 .1379];
den = [1 -1.3746 .6703];

% input and output signals 
N = 1000;
t = [0:N-1]';
unf=randn(N,1)+1;
ynf = dlsim(num ,den,unf);

% The nonlinear element 
znf = (ynf<.5)*.5+((ynf>=.5)&(ynf<=1.5)).*ynf+(ynf>1.5)*1.5;

% 5 Percent process noise and 5 percent measurement noise are added:
u=unf;
usim = unf + .05*std(unf)*randn(size(unf));
y = dlsim(num,den,usim);
zsim = (y<=.5)*.5+((y>.5)&(y<1.5)).*y+(y>=1.5)*1.5;
z = zsim+.05*std(zsim)*randn(size(zsim));

echo off
figure
subplot(1,3,1);plot(100:200,unf(100:200));
axis([100 200 -1 3]);title('noise free input');
subplot(1,3,2);plot(100:200,y(100:200));
axis([100 200 -1 3]);title('lin. outp. incl. proc. noise');
subplot(1,3,3);plot(100:200,z(100:200));
axis([100 200 -1 3]);title('nonlin. outp. incl. noise');
echo on
pause %	Press any key to continue.

clc 
echo off
zd=detrend(z,0);ud=detrend(u,0);
[Sn,R]=dordpi(ud,zd,5);
[Ae,Ce]=destac(R,2);
[Be,De]=destbd(ud,zd,Ae,Ce);
echo on
% Once we know the linear dynamic part, we can estimate the 
% static nonlinear function, between the estimated output of the
% linear model and the output of the wiener system. We can do this
% with for instance an estimation based on Tchebyshev polynomials.
ye=dlsim(Ae,Be,Ce,De,u);
[thl,ze]=chebest(ye,zd,7);

% The next figure shows the relation between the estimated linear
% output y and the nonlinear output z. The x's are the relation
% between the two signals. The line is the approximation with 
% Tchebychev polynomials. The variance accounted for of this 
% estimate is:  

echo off
[thl,ze]=chebest(ye,z,7);
x1=min(ye); x2=max(ye); step=(x2-x1)/N;
x=(x1:step:x2-step)';
line=chebsim(x,thl);
figure
clf
plot(x,line)
hold on
plot(ye,z,'rx')
echo on
vafe=vaf(z,ze)

pause %	Press any key to continue.
clc
% PI-moesp will fail to give a consistent estimate 
% of the linear part unbiased when the input is not Gausian.
% We now need to use a nonlinear optimization technique
% to find the right model. The SMI Toolbox contains a Gauss-Newton
% optimization routine for linear as well as wiener model. It uses the
% Seperable Least Squares technique to minimize the number of
% parameters to be optimized.
% In the next example we have used the same model as before,
% but now a random binary sequence is used as input.
echo off

N=1000;n=2;m=1;l=1;nn=5;
Ts = 0.1;t = [0:Ts:(N-1)*Ts]';
hist1=[];hist2=[];
NUM = [.1578 .1379];
DEN = [1 -1.3746 .6703];

u=0.5+prbn(N,0.2);
usim = u + .05*std(u)*randn(size(u));		% 5 percent process noise 
y = dlsim(NUM,DEN,usim);
zsim = (y<=.5)*.5+((y>.5)&(y<1.5)).*y+(y>=1.5)*1.5;
z = zsim+.05*std(zsim)*randn(size(zsim));	% 5 percent measurement noise

echo on
% The figure shows from left to right, the designed input, the  output
% of the linear part, and the output of the nonlinear part. As can be
% viewed from the plot, the nonlinear element is a saturation at 0.5
% and 1.5. 

echo off
figure
subplot(1,3,1);plot(100:200,u(100:200));
axis([100 200 -1 3]);title('noise free input');
subplot(1,3,2);plot(100:200,y(100:200));
axis([100 200 -1 3]);title('lin. outp. incl. proc. noise');
subplot(1,3,3);plot(100:200,z(100:200));
axis([100 200 -1 3]);title('nonlin. outp. incl. noise');
drawnow
echo on

pause %	Press any key to continue.
clc
% As an initial estimate of the system we use the value PI-moesp
% provides us. Although we know this value will not be
% correct, it is an indication of the true system, and increases
% the chance of finding the global optimum during the
% Gauss-Newton optimization. The gnwisls function only needs the 
% parameters of the linear dynamic part (A,B,C,D) as initial 
% parameters. As with gnlisls, we use d2cm and ss2thon, to get the 
% theta vector.

[Sn,R] = dordpi(u,z,5);
[Ai,Ci] = destac(R,2);
[Bi,Di] = destbd(u,z,Ai,Ci); 

 
% With this initial estimate, we proceed with the Gauss-Newton
% optimization function dslswie. This function takes as input the
% estimated matrices, and the input and output of the system. The next
% input argument is a vector containing some information about the model
% structure.  This vector contains the order nn of Tchebyshev
% polynomials, whether we want to estimate the D matrix and if we need
% an estimate of initial state. In the example we only estimate D. The
% last argument contains options for the optimization routine. dslswie
% uses the function leastsq internally. This function is part of the
% matlab optimization toolbox. 

model=[0 1 0];
options=foptions;
options(1)=1;
%          [A,B,C,D,x0,thl,options] =
%          dslswie(u,z,A,B,C,D,x0,nn,model,options) 

[Ao,Bo,Co,Do,x0o,thlo]=dslswie(u,z,Ai,Bi,Ci,Di,[],7,model,'on',options);

% Finally the variance accounted for is calculated. Let's hope it is
% over 99%.

ye=dlsim(Ao,Bo,Co,Do,u);
ze=tchebsim(ye,thlo);
vaf_end=vaf(z,ze)

echo off
if vaf_end>99
  eval('lalala','')
disp('Hallelujah');
end

echo on
% The result of the optimization can be viewed by ploting the
% nonlinear output as a function of the estimated linear output. 
%We see that the nonlinear element is clearly visible
echo off
figure
x1=min(ye); x2=max(ye); step=(x2-x1)/N;
x=(x1:step:x2-step)';
line=tchebsim(x,thlo);
plot(x,line)
hold on
plot(ye,z,'rx')
axis([-0.5,2,0.2 1.8])
title('nonlinear output vs estimated linear output') 

% This is the end of demo 3.
% Thank you for using SMI.

echo off
% Bert Haverkamp, April 1997
% copyright (c) 1997 Bert Haverkamp










