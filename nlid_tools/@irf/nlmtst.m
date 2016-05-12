function I=nlmtst(i)
% irf/nlmtst
% test IRfs

% one sided first
%
delete(get(0,'children'));
z=nlid_sim ('L1');
figure(1);
i=irf(z,'nLags',101);
plot(i)
figure(2);
r1=nlid_resid(i,z);
title('Residuals');
figure(3);
% Hessian
h=hessian(i, z);
mesh(h);

%
% Convert to fresp;
plot (fresp(i));

%
% Convert to nldat
plot (nldat(abs(i)));

%
% two sided 
%
z1=nldat;
z1(:,1)=z(:,2);
z1(:,2)=z(:,1);
i1=irf(z1,'nSides',2,'nLags',50);
figure(3);
plot(i1)

% Demosntate options for handling noise; 
noise=randvar;
noiseSd=double(std(z(:,2)))/5;
set(noise,'std',noiseSd);
N=nlsim(noise,z(:,2));
zN=z;
zN(:,2)=z(:,2)+N; 


i2=irf(zN,'nLags',101,'irfIdMethod','corr'); 
i3=irf(zN,'nLags',101,'irfIdMethod','pseudo','irfPseudoInvMode','full' );
i4=irf(zN,'nLags',101,'irfIdMethod','pseudo','irfPseudoInvMode','auto' )
figure(99);
i5=irf(zN,'nLags',101,'irfIdMethod','pseudo','irfPseudoInvMode','manual' ,'irfFigNum',99);
figure (2); clf
subplot (2,2,1); plot(i2);title('Correlation Method'); 
subplot (2,2,2); plot(i3);title('Full pseduoinverse'); 
subplot (2,2,3); plot(i4);title('Pseduoinverse - automatic'); 
subplot (2,2,4); plot(i5);title('Pseduoinverse -manual'); 

% 
plot (smo(i,5));




% % test time varying IRF
% % tvtest
 dt=0.01;
 G=randvar;
 nSamp=200;
 nReal=200;
 x=zeros(nSamp,1,nReal);
 X=nldat(x); 
 X=nlsim(G,x);
 set(X,'domainIncr',dt); 
 Y=smo(X,5);
 Z=cat(2,X,Y);

 % TI IRF for each realization 
  iTI=irf(Z,'nSides',2,'nLags',9);
  ypTI=nlsim(iTI,X); 
  v1=vaf(Y,ypTI,'realization');
 plot (v1); 
  
  % Estimate TV IRF for TI case
  
  iTV=irf(Z,'nSides',2,'nLags',9,'tvFlag',true,'irfIdMethod','corr');
  ySimTV = nlsim(iTV, X);
  V= vaf(Y,ySimTV,'realization');
plot(V)  

 %% Generate an TV IRF
 iTV=iTI;
 set(iTV,'tvFlag',true); 
 for i=50:125,
     iTV(:,1,i)=iTI(:,1,i)*.2;
 end
plot(iTV)
yTV=nlsim(iTV,X(:,1,:));

 % TV IRFs
 zTV=cat(2,X,yTV);
 iTVe=irf(zTV,'nSides',2,'nLags',9,'tvFlag',true,'irfIdMethod','pseudo');
 vIRF=vaf(iTV,iTVe,'realization'); 
 plot(vIRF);
 yTVp=nlsim(iTVe,X);
 vPred=vaf(yTV,yTVp,'realization');
plot(vPred) 
   
 %% Problem with scaling in nlsim or NLIDENT?
 
 
 Y
 ;
% y=smo(x,3);
% y(51:100,:)=10*y(51:100,:);
% 
% %
% % Now test nlid tools
% %
% x=reshape (x,100,1,75);
% y=reshape (y,100,1,75);
% z=cat(2,x,y);
% Z=nldat(z);
% 
% i=irf(Z,'nsides',2,'nlags',9,'TV_Flag','Yes','Method','corr');
% ip=zero_pad(i);
% Yp=nlsim(ip,Z(:,1,:));
% figure (4); plot (i);
% title ('Time varying inpulse Response');
