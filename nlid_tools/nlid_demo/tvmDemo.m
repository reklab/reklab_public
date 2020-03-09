function tvmDemo
% tvmDemo - demonstration of tv idetification unsing nlid_tools
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
set(TVirf,'elements',tvI,'tvStart',0,'tvIncr',.0
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
tiIRF=nlident(I,Z(:,:,1));
plot (tiIRF);
figNum=figNum+1; figure(figNum);clf;
[R,V,yp]=nlid_resid(tiIRF,Z(:,:,1));




%% Estimate a TV IRF using ensemble method 
disp('TV IRF Ensemble method');
tvIRFensemble = tvm;
set(tvIRFensemble,'tvIdentMethod','ensemble'); 
tvIRFensemble=nlident(tvIRFensemble, Z, I);
figNum=figNum+1; figure(figNum);clf;

plot(tvIRFensemble);
title('Estimtated TVIRf'); 
figNum=figNum+1; figure(figNum);clf;

title('TV IRF ensemble estimate');
figNum=figNum+1; figure(figNum);clf;

tvResid(tvIRFensemble,Z);


%% Estimate a TV IRF using basis expansion 
nTrials=2; 
disp('Estimate a TV IRF using basis expansion') 
tvIRFbasis = tvm;
set(tvIRFbasis,'tvIdentMethod','basisexpansion'); 
BF=polynom;
set(BF,'polyType','Bspline','polyOrder',10, 'polyRange',[ 0 10],'splineSD',1);
nodeLocations = .5:1:10;
set(BF,'splineCenters',nodeLocations(:));
t=domain(Z);
BF1=basisfunction(BF,t); 
plot(BF1); 
figNum=figNum+1; figure(figNum);clf;
tvIRFbasis=nlident(tvIRFbasis, Z(:,:,1:nTrials), I,BF,'periodic','yes','method','Bayes');

plot(tvIRFbasis);
title(['TV IRF basic function  estimate, Ntrials=' num2str(nTrials)]);
figNum=figNum+1; figure(figNum);clf;
tvResid(tvIRFbasis,Z(:,:,nTrials));
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
set(Zpoly,'chanNames',{ 'Input' 'Output'},'comment', 'TV polynominal');

figNum=figNum+1; figure(figNum);clf;
plot(Zpoly);
streamer('Simulated TV polyonmial response'); 


% Identify a TV polynomial
TVpolyIdent=nlident(tvm,Zpoly,p);
figNum=figNum+1; figure(figNum);clf;
plot(TVpolyIdent);
title('Estimated TV Polynomial'); 
%
simY=nlsim(TVpolyIdent,X);
figNum=figNum+1; figure(figNum);clf;
subplot(3,1,1);
plot(Ypoly);
title('Output');
subplot (3,1,2);
plot (simY);
title('Estimated Output');subplot(3,1,3);
plot(Ypoly-simY);
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
plot(Znlbl)


TVnlblID=nlident(tvm,Zpoly,NLBL);
yPre=nlsim(TVnlblIdent, X);












   