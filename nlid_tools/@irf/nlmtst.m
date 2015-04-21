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

%
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
% dt=0.01;
% x=randn(100,75);
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
% set(Z,'domainincr',.01);
% i=irf(Z,'nsides',2,'nlags',9,'TV_Flag','Yes','Method','corr');
% ip=zero_pad(i);
% Yp=nlsim(ip,Z(:,1,:));
% figure (4); plot (i);
% title ('Time varying inpulse Response');
