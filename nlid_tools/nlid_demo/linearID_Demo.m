function irfDemo (inputCutOff, noiseLevel, nLags,nFFT,  nSamp, sysType)
% mod9Demo (inputCutOff, noiseLevel, nLags, nSamp)
% Demonstrate IRF estimation for various systems
% inputCutOff - normalized cutoff for input singal (0-1)
% noiseLevel - ration of noise STD to output SRD (0 -100);
% nLags - lenght of IRF
% nSamp - number of samples.
% sysType - system to simulate

%%  IRF demos
if nargin < 1,
    inputCutOff=1;
end
if nargin <2
    noiseLevel=0;
end
if nargin <3
    nLags=50;
end
if nargin <4
    nFFT=100;
end
if nargin <5
    nSamp=20000;
end
if nargin <6
    sysType='all';
end

delete(get(0,'children'));
set(0,'DefaultFigureWindowStyle','docked')

nSides=2;

%% Generate input signal
u=randn(nSamp,1);
if inputCutOff<1,
    [b,a]=butter(2,inputCutOff/2, 'low');
    u=filter(b,a,u);
end
figure(1);
U=nldat(u,'domainIncr',.01,'comment','Input');
subplot (2,2,1);
plot (U(1:1000));
subplot (2,2,2);
SU=spect(U,'nFFT',100);
set(U,'comment','Input Spectrum','domainName','sec');
plot (SU);title('Spectrum');
fig_mod(1,'title_size',14);
subplot (2,2,2);
CU=cor (U,'nSides',2);
plot(CU);title('Autocorrelation');
subplot (2,2,[3 4]);
plotHessian(U);title('Hessian');
streamer(['Input signal BW=' num2str(inputCutOff)],.9);

%
% Generare noise signal
rNoise = randvar;
noise =nlsim(rNoise,domain(U));
noise=noise;
switch sysType
    case {'static' , 'all'}
        %% Static Linear
        [z,m]=nlid_sim('static_linear',U,'noise_level',0);
        nSides=2;
        titleStr='static Linear' ;
        plotIRF; plotFR
        %% Dynamic LowPass
    case {'lowpass' , 'all'}
        z=nlid_sim('L1',U,'noise_level', 0);
        nSides=1;
        titleStr='Dynamic LowPass' ;
        plotIRF; plotFR
        %% Dynamic HighPass
    case {'highpass' ,'all'}
        z=nlid_sim('H1',U,'noise_level', 0);
        nSides=2;
        titleStr='Dynamic high Pass' ;
        plotIRF; plotFR
        %% Static Linear with Delay
    case {'static delay' , 'all'}
        
        z=nlid_sim('static_linear',U,'noise_level', 0,'delay_time', .100);
        nSides=1;
        titleStr='Statc Linear with Delay' ;
        plotIRF; plotFR
        %% LowPass with Delay
    case {'lowpass_delay' , 'all'}
        z=nlid_sim('L1',U,'noise_level', 0,'delay_time', .100);
        nSides=1;
        titleStr='Low Pass with Delay' ;
        plotIRF; plotFR
        %%
    case {'Hammerstein' , 'all'}
        z=nlid_sim('N3L',U,'noise_level',0);
        nSides=1;
        titleStr='Hammerstein System ' ;
        plotIRF; plotFR
    otherwise
        error( 'specified sysType do not exist');
end
    function plotIRF
        figure(3);
        disp(titleStr);
        fullTitle=[titleStr '; Input BW= ' num2str(inputCutOff) '; Noise:' num2str(noiseLevel) '; nSamp=' num2str(nSamp)];
        set (gcf,'name',fullTitle);
        zIn=z(:,1);acIn=cor(zIn,'nSides',2,'nLags',nLags);
        % Add noise
        stdZout=std(double(z(:,2)));
        stdNoise=std(double(noise));
        gain =noiseLevel*stdZout/stdNoise;
        zOut=z(:,2)+(noise*gain);
        Z=cat(2,zIn,zOut);
        acOut=cor(zOut,'nSides',2,'nLags',nLags);
        xCor=cor(Z,'nSides',nSides,'nLags',nLags);
        I=irf(Z,'nSides',nSides,'nLags',nLags);
           Ia=irf(z,'nSides',nSides,'nLags',nLags,'irfPseudoInvMode','manual');
        zPre=nlsim(I,zIn);
        subplot (2,3,2);
        plot (acIn);
        title('Input autocorrelation');
        
        subplot (2,3,3);
        plot (acOut);
        title('Output autocorrelation');
        subplot (2,3,4);
        plot (xCor);
        title('Cross Correlation');
        subplot (2,3,5);
        plot (I)
        subplot (2,3,6);
     
        plot (Ia);
        title('pseudoInverse IRF');
        figure(4)
        
        iVAF=double(vaf(zOut,zPre));
        % Plot input, output, predicted output
        zPlot = cat(2,Z, zPre);
        subplot (4,1,1);
        plot (zPlot(50:550,1));
        title (['VAF=' num2str(iVAF)]);
        ylabel('Input');
        xlabel('');
        
        subplot (4,1,2);
        plot (zPlot(50:550,2),'line_color','r');
        ylabel('Output');
        title('');
        xlabel('');
        
        subplot (4,1,3);
        plot (zPlot(50:550,3),'line_color','g');
        ylabel('Predicted');
        title('');
        subplot (4,1,4);
        plot (zPlot(50:550,2)-zPlot(50:550,3),'line_color','g');
        ylabel('Residuals');
        title('');
          streamer(fullTitle,.90);
        fig_mod(4,'title_size',12);
        
    end
function plotFR
z1=z;
figure(5);
S=spect(z1,'nFFT',nFFT);
plot (abs(S),'xmode','log','ymode','db');
title('Input Spectrum'); 
set (gcf,'name','Input Spectrum');

figure(6);clf;
fullTitle=[titleStr '; Input BW= ' num2str(inputCutOff) '; Noise Level:' num2str(noiseLevel)];
set (gcf,'name',fullTitle);
set(z1,'comment',fullTitle);
F=fresp(z1,'nFFT',nFFT);
plot(F);
title(fullTitle);
end
end
function plotHessian(U)
c=cor(U,'nLags',32);
T=toeplitz(double(c));
H=T'*T;
mesh(H)
condNum=cond(H);
title(['Hessian. Condition number=' num2str(condNum)]);
end
