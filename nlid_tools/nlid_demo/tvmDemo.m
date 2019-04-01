function tvmDemo
% tvmDemo - demonstrate tvm propoerties
%% TV IRF
%%Generate an TV IRF
figNum=0;
disp ('tvIRF identification');
disp('Generate TV IRF');
gain=ones(200,1);
gain=cat(1,gain,linspace(0,5,600)');
gain=cat(1,gain, ones(200,1));
plot(gain)
tvI={};
 for i=1:1000,
     tvTemp=irf2(irf,'g',gain(i));
     tvI{i}=tvTemp;
 end
 TVirf=tvm;
set(TVirf,'elements',tvI,'tvStart',0,'tvIncr',.01);
figNum=figNum+1; figure(figNum);clf;
plot(TVirf); title('Simuulated TV IRF'); 
%% simulate TV responses
disp('Simulate at TV Response');
x=randn(1000,1,200);
X=nldat(x,'domainIncr',.01);
Y=nlsim(TVirf,X);
Z=cat(2,X,Y); 
set(Z,'chanNames',{ 'Input' 'OutPut'},'comment','Simulated TV trial');
figNum=figNum+1; figure(figNum);clf;

plot(Z);

%% Estimate a TI IRF
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
tvIRFensemble = tvm;
set(tvIRFensemble,'tvIdentMethod','ensemble'); 
tvIRFensemble=nlident(tvIRFensemble, Z, I);
figNum=figNum+1; figure(figNum);clf;

plot(tvIRFensemble);
figNum=figNum+1; figure(figNum);clf;

title('TV IRF ensemble estimate');
figNum=figNum+1; figure(figNum);clf;

tvResid(tvIRFensemble,Z);


%% Estimate a TV IRF using basis expansion 
disp('Estimate a TV IRF using basis expansion') 
tvIRFbasis = tvm;
set(tvIRFbasis,'tvIdentMethod','basisexpansion'); 
BF=polynom;
set(BF,'polyType','B_spline','polyOrder',10, 'polyRange',[ 0 10],'B_spline_SD',1);
nodeLocations = .5:1:10;
set(BF,'polyCoef',nodeLocations(:));
t=domain(Z);
BF1=basisfunction(BF,t); 
plot(BF1); 
figNum=figNum+1; figure(figNum);clf;
tvIRFbasis=nlident(tvIRFbasis, Z(:,:,1:50), I,BF,'periodic','yes','method','Bayes');
figNum=figNum+1; figure(figNum);clf;

plot(tvIRFbasis);
title('TV IRF basic function  estimate');
figNum=figNum+1; figure(figNum);clf;
tvResid(tvIRFbasis,Z);
subplot (3,1,1); title('Residuals for TV IRF Basis Function Estimate');


%% TV Polynomial
disp('TV polynomial demo');
% Generate a TV Polynomial
p=polynom('polyType','power','polyOrder',2,'polyOrderMax',3)
PI=[];
coef=[ 0 1 1];
for i=1:1000,
    coef(2)=1+gain(i);
     pTemp=set(p,'polyCoef',coef);
     PI{i}=pTemp;
end
 TVpoly=tvm;
 set(TVpoly,'elements',PI,'tvIncr',X.domainIncr);
figNum=figNum+1; figure(figNum);clf;

 plot(TVpoly);
 title('TV Polyomial'); 
 % Simulate response to TV polynomial
Ypoly=nlsim(TVpoly,X);
Zpoly=cat(2,X,Ypoly); 
figNum=figNum+1; figure(figNum);clf;
plot(Zpoly);
title('Simulated TV polynmial response'); 


% Identify a TV polynomial
TVpolyIdent=nlident(tvm,Zpoly,p);
figNum=figNum+1; figure(figNum);clf;
plot(TVpolyIdent);
title('Estimated TV Polynomial'); 
%
simY=nlsim(TVpolyIdent,X);
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
%% TV Hammerstein Identification - ensemble method
disp('TV Hammerstein  - In progres');
polyElements=TVpoly.elements;
irfElements=TVirf.elements;
nElements=length(polyElements);
NLBL=nlbl;

for i=1:1000,
    set(NLBL,'elements',{polyElements{i} irfElements{i}});
    nlblElements{i}=NLBL;
end
TVnlbl=tvm;
set(TVnlbl,'elements',nlblElements,'tvStart',0,'tvIncr',.01);
y1=nlsim(TVirf,Ypoly);

yNLBL=nlsim(TVnlbl,X);

Znlbl=cat(2,X,yNLBL); 
figNum=figNum+1; figure(figNum);clf;
plot(simNLBL)


TVnlblIdent=nlident(tvm,Zpoly,NLBL);
yPre=nlsim(TVnlblIdent, X);












   