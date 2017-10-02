function tvmDemo
% tvmDemo - demonstrate tvm propoerties
%% TV IRF
%%Generate an TV IRF
figNum=0;
disp ('tvIRF identification');
disp('Generate TV IRF');
gain=zeros(200,1);
gain=cat(1,gain,linspace(0,5,600)');
gain=cat(1,gain, zeros(200,1));
plot(gain)
tvI={};
 for i=1:1000,
     tvTemp=irf2(irf,'g',gain(i));
     tvI{i}=tvTemp;
 end
 TVM=tvm;
set(TVM,'elements',tvI,'tvStart',0,'tvIncr',.01);
figNum=figNum+1; figure(figNum);clf;
plot(TVM);
%% simulate TV responses
disp('Simulate at TV Response');
x=randn(1000,1,200);
X=nldat(x,'domainIncr',.01);
Y=nlsim(TVM,X);
Z=cat(2,X,Y); 
set(Z,'chanNames',{ 'Input' 'OutPut'});
figNum=figNum+1; figure(figNum);clf;

plot(Z);

%% Estimate at TI IRF
disp('TI IRF');
figNum=figNum+1; figure(figNum);clf;
I=irf;
set(I,'nLags',100);
tiIRF=nlident(I,Z);
[R,V,yp]=nlid_resid(tiIRF,Z);
figNum=figNum+1; figure(figNum);clf;
subplot (2,1,1);
plot(V);
title('TI IRF VAF');
xlabel('Realization');
subplot (2,1,2);
plot(R);
title('TI Residuals');
xlabel('Time (s)');






%% Estimate a TV IRF using ensemble method 
disp('TV IRF Ensemble method');
tvIRF = tvm;
set(tvIRF,'tvIdentMethod','ensemble'); 
tvIRF=nlident(tvIRF, Z, I);
figNum=figNum+1; figure(figNum);clf;

plot(tvIRF);
figNum=figNum+1; figure(figNum);clf;

title('TV IRF ensemble estimate');
figNum=figNum+1; figure(figNum);clf;

tvResid(tvIRF,Z);


%% Estimate a TV IRF using basis expansion 
disp('Estimate a TV IRF using basis expansion') 
tvIRF = tvm;
set(tvIRF,'tvIdentMethod','basisexpansion'); 
BF=polynom;
set(BF,'polyType','B_spline','polyOrder',10, 'polyRange',[ 0 10],'B_spline_SD',1);
nodeLocations = .5:1:10;
set(BF,'polyCoef',nodeLocations(:));
t=domain(Z);
BF1=basisfunction(BF,t); 
plot(BF1); 
figNum=figNum+1; figure(figNum);clf;
tvIRF=nlident(tvIRF, Z(:,:,1:50), I,BF,'periodic','yes','method','Bayes');
figNum=figNum+1; figure(figNum);clf;

plot(tvIRF);
title('TV IRF basic function  estimate');
figNum=figNum+1; figure(figNum);clf;
tvResid(tvIRF,Z);
subplot (3,1,1); title('Residuals for TV IRF Basis Function Estimate');


%% TV Polynomial
disp('TV polynomial demo');
% Generate a TV PolynomiAL
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
figNum=figNum+1; figure(figNum);clf;

 plot(TVM);
 title('TV Polyomial'); 
 % Simulate response to TV polynomial
Y=nlsim(TVM,X);
Z=cat(2,X,Y); 
figNum=figNum+1; figure(figNum);clf;
plot(Z);
title('Simulated TV polynmial response'); 


%% Identify a TV polynomial
tvIdent=nlident(tvm,Z,p);
figNum=figNum+1; figure(figNum);clf;
plot(tvIdent);
title('Estimated TV Polynomial'); 
%
simY=nlsim(tvIdent,X);
figNum=figNum+1; figure(figNum);clf;
subplot(3,1,1);
plot(Y);
title('Output');
subplot (3,1,2);
plot (simY);
title('Estimated Output');
subplot(3,1,3);
plot(Y-simY);
title('Residuals');
%% TV Hammerstein Identification 
disp('TV Hammerstein  - Still to come');




   