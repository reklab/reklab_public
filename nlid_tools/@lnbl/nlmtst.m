function i=lnmtst(i)
% lnlbl/nlbl
% test of lnbl identification
%
clear all
disp('nlmtst for lnbl');

z=nlid_sim ('LN3');
i=lnbl;
%
% Generate a ln system;

I=irf2(irf);
P=polynom('polyType','power','polyCoef', [ 100 10 1],'polyRange',[-10 10]);
plot(P)
W = lnbl('elements', { I  P});
plot(W) 
X=10*nldat(randn(1000,1),'domainIncr',.01);
Y=nlsim(W,X);
Z=cat(2,X,Y);
i=nlident(lnbl,Z,'lnIdMethod','bussgang','nLags',128);
figure(1);plot(i);
i=nlident(lnbl,Z,'lnIdMethod','hk','nLags',128);
figure(2);plot(i);
i=nlident(lnbl,Z,'lnIdMethod','phk','nLags',128);
figure(3);plot(i);
i=nlident(lnbl,Z,'lnIdMethod','lm','nLags',128);
figure(2);plot(i);



x=z(:,1);
y=z(:,2);
yp=nlsim(i,x);
z(:,1)=yp;
figure(2);
plot(z);
vaf(y,yp);

% lnbl/nlmtst
