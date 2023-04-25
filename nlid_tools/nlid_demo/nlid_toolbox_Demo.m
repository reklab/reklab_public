%%nlid_toolbox_demo
% Object oriented toolbox for nonparamtric identification of linear and
% nonlinear systems.
%% Enable cell mode and go through one cell at a time
%$ Getting Help 

%% Display a list of all the classess defined in the nlid_toolbox' 
nlid_help ('classList');
%% Show Properties defined for a class
properties(nldat)
%% show all defined methods for a class nldat
methods(nldat)
%% Show help for one method for  for a class
help nldat/plot
%% Show help summary for all methods for class nldat
disp('Display help on all methods for a class') 
nlid_help('nldat');
%% nldat class - data classcd('
% Create a nldat vector
x=randn(5000,1);
X=nldat(x); % Convert real to nldat
y=double(X); % Convert nldat to real. 

%% 
disp('All properties are defined by default')
disp(X)
%% 
disp('Recall one property Value using "." format')
incr=X.domainIncr;
%% 
disp('Alternatively use the get command')
incr=get(X,'domainIncr');
disp(incr);
%% Set one property value using "." format
X.domainIncr=.001;
disp(X);

%% Set multiple properties using the set command
set(X,'domainIncr',.01,'chanNames',{'X'}, 'comment','Test signal for demonstration purposes', ...
 'chanUnits',{'V'} ); % set domain increment
disp(X) % 
clf
plot(X)

%%nldat object can have multiple channels

Y =smo(X,5) % Smooth X five times
set(Y,'chanNames',{'Y'});
Z=cat(2,X,Y);
set(Z,'comment','nldat with two channels');
clf
plot(Z);
disp(Z)

%%nldat object can have multiple realizatios
Z2=2*Z+5;
Z3=cat(3,Z,Z2);
Z3.comment='Demo signal with 2 channels and 2 realizations'; 
plot(Z3,'realizationMode','offset')
disp(Z3)
% Show plot options
help nldat/plot

%% Use matlab idexing to return samples, channels or realizaions;
%% nldat( iSample, iChannel, iRealization)
 clf; plot(Z3(1:10,2,1))  % Ten samples from one channel and one realization
 clf; plot(Z3(:, 1,1)); title('One channel') % All samples from one channel 
clf;plot( Z3(:,1,:))  % All realizations of one channell
 

%% Aritmetic operations on nldat objects follows  matlab syntax 

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
plot (X); title('X');
subplot (3,1,2);
plot (Y);title('Y');
subplot (3,1,3);title('Z=X+Y');
plot (Z);

%% Many matlab operators are overlaid for nldat objects
methods(nldat);
mean(Z)
mean(Z2)
mean(Z3)%%


%%   pdf - probabilty distribution function class for NLID toolbox
% Compute and plot probability distribution from nldat object
pX = pdf(X);
clf; plot(pdf(X));
pX=pdf(X,'nBins',20,'pdfType','Frequency');
clf; plot (pX);
%% randvar - nlid class supporting random variates
r=randvar;
disp(r)
% See what types are available
set(r,'randvarType');
set(r,'randvarType','Normal');disp(r);
% pdf supports randvar objects as well. 
clf;plot (pdf(r));

%% Generate a realization of 1000 points
% nlsim method of ranvar generates a realization of a randvar object
X=nlsim(r,[1:1000]')
set(X,'domainIncr',.001);
plot (X)
% compare theoretical and observed PDFs
plot (pdf(X,'nBins',20))
px=pdf(r);  % Compute theoretical distribution
h=line(px);set(h,'color','r','linewidth',2)
legend ('data','theoretical');

%% spect - power spectrum class for NLID toolbox.  
% Compute and plot spectrum
X=nlsim(r,[1:10000]')
set(X,'domainIncr',.01);
fX=spect(X,'nFFT', 100);
disp(fX);
subplot (2,1,1);title('Spectrum of X: linear plot'); 
plot (fX);
subplot (2,1,2);
plot(fX,'xmode','log','ymode','db'); title ('Spectrum of X: log-DB plot'); 
%% Compute and plot autocorrelation 
% Default parameters
C=cor(X);
subplot (3,1,1);
plot(C); title('Autocorrelation of X: Default Paramters');
% Nsides =2 ;
subplot(3,1,2);
plot(cor(X,'nSides',2));title('Autocorrelation of X: 2 Sided');
% Nsides =2 nlags=50;;
subplot(3,1,3);
plot(cor(X,'nSides',2,'nLags',50)); title('Autocorrelation of X: 2 Side 100 lags');

%% Simulate a second order low pass system
Z = nlid_sim('L1',X);
set(Z,'comment','Simulated low-pass system'); 
plot(Z); 

%% Correlation properties
% Note that nldat objects are indexed in the same way as matlab variables
% Input Autocorrelation
nLags=30;
subplot(3,1,1);
cX = cor(Z(:,1),'nSides',2,'nLags',nLags,'comment','Input autocorrelation');
plot(cX); 
% output autocorrelation
subplot (3,1,2);
cY = cor(Z(:,2),'nLags',nLags, 'nSides',2, 'comment','Output autocorrelation');
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
plot(sX); title('Input Spectrum'); 
% output spectrum
subplot (3,1,2);
sY = spect(Z(:,2),'nFFT',nFFT, 'comment','Output spectrum');
plot(sY);title('Output Spectrum')
% Input/output cross correlation
subplot(3,1,3);
sXY = spect (Z,'nFFT',nFFT,'comment','Input-output cross-sprectrum'); 
plot (abs(sXY));title('Magnitude of input-output cross spectrum'); 

%% Impulse Response function
plot(Z); 
I=irf(Z,'nLags',51);
clf; plot(I);
%% Simulate response
xIn=Z(:,1);
yOut=Z(:,2);
yPre = nlsim(I,xIn);
plot(yPre);title('Predicted output') 
yResid = yOut - yPre;
set(yResid,'chanNames','Residual'); 
V=vaf(yOut,yPre); 
disp(['Variance accounted for:' num2str(double(V))])
Zp=cat(2,yOut, yPre);
Zp=cat(2,Zp, yResid);
set(Zp,'comment','Resdiual Analysis');
plot(Zp,'nv',3); title('Residuals'); 
streamer(['Variance accounted for:' num2str(double(V))])


%% Frequency Response
F=fresp(Z, 'nFFT',50);
plot(F)
figMod(1,'lineWidth',2);
%% Predict response
yPre = nlsim(F,xIn);
plot(yPre);
yResid = yOut - yPre;
set(yResid,'chanNames','Residual'); 

V=vaf(yOut,yPre); 
disp(['Variance accounted for:' num2str(double(V))])
Zp=cat(2,yOut, yPre);
Zp=cat(2,Zp, yResid);
Zp.comment='Frequency Response Prediction';
plot(Zp,'nv',3)

%% Convert IRf to frequency response
fI = fresp(I);
fI.comment='Frequecy Response from IRF';
plot(fI);

%% Convert  frequency response to IRF 
iFreq = irf(F);
iFreq.comment='IRF from Frequency Response'; 
clf;plot(iFreq);







