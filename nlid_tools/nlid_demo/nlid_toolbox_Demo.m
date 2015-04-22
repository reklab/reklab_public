%%nlid_toolbox_demo
%
% Script to demonstrate basic functionality of nlid_tools
%
%
% Enable cell mode and go through one cell at a time
%$ Getting Help 

%% Display a list of all the classess defined in the nlid_toolbox' 
nlid_help ('classList');
%% show all defined methods for class nldat
methods(nldat)
%% Show help for one methods for a class
help nldat/plot
%% Show help summary for all methods for a class nldat
disp('Display help on all methods for a class') 
nlid_help('nldat');
%% nldat class - data class
x=randn(5000,1);
X=nldat(x); % Convert real to nldat
% set(X,'domainIncr',.01,'chanNames',{'X'}, 'comment','Test signal for demonstration purposes', ...
% 'chanUnits',{'V'} ); % set domain increment
disp(X) % Show properties
clf
figure(1); clf
set(X,'domainIncr',.000001)

plot(X)

%% Arithmetic
clf;
Y= 2*X;
set (Y,'comment', ' 2 * X','chanNames',{'Y'});
subplot (2,1,1);
plot (X);
subplot (2,1,2);
plot (Y)
%% Summation
Z = X+Y;set(Z,'chanNames',{'Z'},'comment',{'X+Y+Z'});
subplot (3,1,1);
plot (X);
subplot (3,1,2);
plot (Y);
subplot (3,1,3);
plot (Z);

%% Array operations
Z=cat(2,X,Y);
plot (Z)
%% indexing
clf; plot (Z(1:100,:));

%%
clf; plot (Z(1:100,1))


%% Compute and plot probability distribution
pX = pdf(X);
plot(pdf(X));
pX=pdf(X,'nBins',20,'pdfType','frequency');
plot (pX);
%% randvar
r=randvar;
disp(r)
% See what types are available
set(r,'randvarType');
set(r,'randvarType','Uniform');disp(r);
clf;plot (pdf(r));

% Generate a realization of 1000 points
x=nlsim(r,[1:1000]')
set(x,'domainIncr',.001);
plot (x)
% compare theoretical and observed PDFs
plot (pdf(x,'nBins',20))
px=pdf(r);
h=line(px);set(h,'color','r','linewidth',2)


%% Compute and plot spectrum
fX=spect(X);
disp(fX);
subplot (2,1,1);
plot (fX);
subplot (2,1,2);
plot(fX,'xmode','log','ymode','db')
%% Compute and plot autospectrum 
% Default parameters
C=cor(X);
subplot (3,1,1);
plot(C);
% Nsides =2 ;
subplot(3,1,2);
plot(cor(X,'nSides',2));
% Nsides =2 nlags=50;;
subplot(3,1,3);
plot(cor(X,'nSides',2,'nLags',50));

%% Simulate a second order low pass system
Z = nlid_sim('L1',X);
plot(Z)

%% Correlation properties
% Note that nldat objects are indexed in the same way as matlab variables
% Input Autocorrelation
nLags=30;
subplot(3,1,1);
cX = cor(Z(:,1),'nLags',nLags,'comment','Input autocorrelation');
plot(cX); 
% output autocorrelation
subplot (3,1,2);
cY = cor(Z(:,2),'nLags',nLags, 'comment','Output autocorrelation');
plot(cY);
% Input/output cross correlation
subplot(3,1,3);
cXY = cor (Z,'nLags',nLags,'comment','Input-output cross-correlation'); 
plot (cXY);

%% Spectral  properties
% 
% Input spectrum 
nFFT=200;
subplot(3,1,1);
sX = spect(Z(:,1),'nFFT',nFFT,'comment','Input Spectrum');
plot(sX); 
% output specrum
subplot (3,1,2);
sY = spect(Z(:,2),'nFFT',nFFT, 'comment','Output spectrum');
plot(sY);
% Input/output cross correlation
subplot(3,1,3);
sXY = spect (Z,'nFFT',nFFT,'comment','Input-output cross-sprectrum'); 
plot (sXY);

%% Impulse Response function
I=irf(Z,'nLags',51);
plot(I);
%% Simulate response
xIn=Z(:,1);
yOut=Z(:,2);
yPre = nlsim(I,xIn);
plot(yPre);
yResid = yOut - yPre;
set(yResid,'chanNames','Residual'); 

V=vaf(yOut,yPre); 
disp(['Variance accounted for:' num2str(double(V))])
Zp=cat(2,yOut, yPre);
Zp=cat(2,Zp, yResid);
plot(Zp,'nv',3)
streamer(['Variance accounted for:' num2str(double(V))])


%% Frequency Response
F=fresp(Z, 'nFFT',50);
plot(F)
%% Predict response
yPre = nlsim(F,xIn);
plot(yPre);
yResid = yOut - yPre;
set(yResid,'chanNames','Residual'); 

V=vaf(yOut,yPre); 
disp(['Variance accounted for:' num2str(double(V))])
Zp=cat(2,yOut, yPre);
Zp=cat(2,Zp, yResid);
plot(Zp,'nv',3)

%% Convert IRf to frequency response
fI = fresp(I);
plot(fI);

%% Convert  frequency response to IRF 
iFreq = irf(F);
plot(iFreq);







