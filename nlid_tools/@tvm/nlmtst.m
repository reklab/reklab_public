function [MP, MI, Z, yp] = nlmtst (TVM)
% test function for tvm

% Copyright 2000, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

figure(1)
x= linspace (-3,3);
xr=reshape(repmat(x,40,1),40,1,100);
y1=10+20*x + 3*x.^2;
y1r=reshape(repmat(y1,20,1),20,1,100);
y2=5 -10*x +5*x.^3;
y2r=reshape(repmat(y2,20,1),20,1,100);
y=cat(1,y1r,y2r);
z=cat(2,xr,y);
Z=nldat(z);
figure(1);
plot(Z,'moder','mesh');
set(Z,'domainincr',.01,'domainstart',0,'domainname','Seconds');
M=tvm;
set(M,'Model_Type','polynom');
M=nlident(M,Z);
yp=nlsim(M,Z(:,1,:));
disp('polynom done')
figure(2);
plot (M)
MP=M;
%
% Test for tv irf
% test time varying IRF
% tvtest

dt=0.01;
x=randn(100,75);
y=smo(x,3);
y(51:100,:)=10*y(51:100,:);

x=reshape (x,100,1,75);
y=reshape (y,100,1,75);
z=cat(2,x,y);
Z=nldat(z);
set(Z,'domainincr',.01);
M=tvm;
set (M,'Model_Type','irf');
set(M,'nsides',2,'nlags',9,'Method','pseudo');
MI=nlident(M,Z)
Yp=nlsim(MI,Z(:,1,:));
figure (4); plot (MI);
title ('Time varying inpulse Response');
