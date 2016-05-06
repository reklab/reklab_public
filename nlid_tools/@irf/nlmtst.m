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
%
% two sided 
%
z1=nldat;
z1(:,1)=z(:,2);
z1(:,2)=z(:,1);
i1=irf(z1,'nSides',2,'nLags',50);
figure(3);
plot(i1)

% Test options
noise=randvar;
set(noise,'std',double(std(z(:,2))));

i2=irf(z,'nLags',101,'irfIdMethod','corr');
i3=irf(z,'nLags',101,'irfIdMethod','pseudo','irfPseudoInvMode','full' );
i4=irf(z,'nLags',101,'irfIdMethod','pseudo','irfPseudoInvMode','auto' );
i5=irf(z,'nLags',101,'irfIdMethod','pseudo','irfPseudoInvMode','manual' );



% 
plot (smo(i,5));


% Hessian
h=hessian(i, z);

%
% Convert to fresp;
plot (fresp(i));

%
% Convert to nldat
plot (nldat(abs(i)));



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
  ypTI=nlsim(iTI,Z(:,1,:));
  vaf(Y,ypTI)
  
 % Generate TV IRF
 iTV=iTI;
 set(iTV,'tvFlag',true); 
 for i=50:125,
     iTV(:,1,i)=iTI(:,1,i)*0;
 end
yTV=nlsim(iTV,X(:,1,1));

 % TV IRFSS
   iTV=irf(Z,'nSides',2,'nLags',9,'tvFlag','Yes','irfIdMethod','corr');
 
 
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
