function tvmDemo
% tvmDemo - demonstrate tvm propoerties
%% TV IRF
%Generate an TV IRF
gain=ones(200,1);
gain=cat(1,gain,linspace(1,5,600)');
gain=cat(1,gain, ones(200,1));
plot(gain)
tvI={};
 for i=1:1000,
     tvTemp=irf2(irf,'g',gain(i));
     tvI{i}=tvTemp;
 end
 TVM=tvm;
 
set(TVM,'elements',tvI,'tvStart',0,'tvIncr',.01);
figure(1);clf;
plot(TVM);

x=randn(1000,1,200);
X=nldat(x,'domainIncr',.01);

Y=nlsim(TVM,X);
Z=cat(2,X,Y); 

I=irf;
set(I,'nLags',50);
tvIdent=nlident(tvm,Z,I);

simY=nlsim(tvIdent,X);

figure(2);clf
subplot(3,1,1);
plot(Y);
title('Output');
subplot (3,1,2);
plot (simY);
title('Estimated Output');
subplot(3,1,3);
plot(Y-simY);
title('Residuals');


%% TV Polynomial

p=polynom('polyType','power','polyOrder',2)
PI=[];
coef=[ 0 1 1];
for i=1:1000,
    coef(3)=gain(i);
     pTemp=set(p,'polyCoef',coef);
     PI{i}=pTemp;
end
 TVM=tvm;
 set(TVM,'elements',PI,'tvIncr',X.domainIncr);
 figure (3);clf
 plot(TVM);
 title('TC Polyomial'); 
 
Y=nlsim(TVM,X);
Z=cat(2,X,Y); 


tvIdent=nlident(tvm,Z,p);

simY=nlsim(tvIdent,X);

figure(4);clf
subplot(3,1,1);
plot(Y);
title('Output');
subplot (3,1,2);
plot (simY);
title('Estimated Output');
subplot(3,1,3);
plot(Y-simY);
title('Residuals');




; 
   