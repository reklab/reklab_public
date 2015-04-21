% Demo_502  - demosntrationfor 502

%% Tool box Demo
nlid_toolbox_Demo
%% Linear identification 
% White Input , no noise

inputBW  = 1.;
noiseLevel = 0;
irfLen = 50;
nFFT = 100;
dataLen =1000;
sysNum=2;

% Linear system identification 
while 1
    
sysList= { 'static' 'lowpass' 'highpass' 'static delay' 'lowpass delay' 'Hammerstein' 'done'}
sysNum=menu('Select system',sysList);
sysType=sysList{sysNum};
if strcmp('sysType','done'),
    break;
end

inputBW = input_d ('Input bandwith  (0-1)', inputBW,0, 1);
noiseLevel= input_d('Noise level (0 - 100)', noiseLevel,0,1000);
irfLen=input_d('IRF length (0-1000)', irfLen,10,1000);
nFFT = input_d('FFT Length',nFFT, 10,10000);
dataLen = input_d('Data length' ,dataLen, 100,100000);
linearID_Demo ( inputBW, noiseLevel,  irfLen, nFFT, dataLen, sysType);
end


%
%% nonlinear identification
%
r=randvar;
x=nlsim(r,[0:10E4]');set(x,'domainIncr',.01);
noiseLevel=0;

% Weiner system 
nonlinearID_demo ('LN3', 10*x, noiseLevel);
% Hammertstein System
nonlinearID_demo ('N3L', 10*x, noiseLevel);
%
% Now add some noise
%

noiseLevel=1;
nonlinearID_demo ('LN3', 10*x, noiseLevel);

% Band limit the input
xBand=smo(x,2); noiseLevel=0;
xRange=range(10*x);
xBandRange=range(xBand);
xBand = xBand*double(xRange)./double(xBandRange);

nonlinearID_demo ('LN3', xBand, noiseLevel);
nonlinearID_demo ('N3L', xBand, noiseLevel);



 

