function noise = snr_tv_gen(y,SNR,wL,nType,Fs)

%% Ehsan Sobhani Tehrani, Sept. 9, 2015
%++ This function generates a non-stationary noise signal that generates a fixed SNR for a the non-stationary output
%-- y: The noise-free signal
%-- SNR: Signal to Noise Ratio
%-- wL: Window Length in samples
%-- nType: The type of noise
%-- Fs: Sampling Time ('white','colored',experimental')

%-- noise: The generated non-stationary noise

%% List of Modifications:
%-- I added an option to generate white or colored noise, nType
%-- I added an input for the sampling time of the noise (needed for filter design in the case of colored noise)
%-- Jan. 18, 2016: I added the option for generating a TV experimental noise  

N = length(y);   %-- number of samples

noise = zeros(N,1);

nW = fix(N/wL);

switch nType
    case 'white'
        for i = 1:nW
            indx1 = (i-1)*wL + 1;
            indx2 = indx1 + wL - 1;
            noise(indx1:indx2,1) = snr_gen(y(indx1:indx2,1),SNR);
        end            
    case 'colored'
        nFpass = 15; %12; %12; %15;  %-- Hz
        nFstop = 20; %17; %18; %20;  %-- Hz
        Ap = 1; %0.2;
        Ast = 60; %40;
        dFilter = fdesign.lowpass('Fp,Fst,Ap,Ast',nFpass,nFstop,Ap,Ast,Fs);
        hFilter = design(dFilter,'butter');
        % hFilter = design(dFilter,'equiripple','StopbandShape','linear','StopbandDecay',20);
        
        % filtOrder = 6;
        % dFilter = fdesign.lowpass('N,Fp,Ap,Ast',filtOrder,nFpass,Ap,Ast,Fs);
        % hFilter = design(dFilter);      
        
        nWhite = 1*randn(N,1);
        nColored = filter(hFilter,nWhite); %-- Making the noise colored
        %-- Scaling the colored noise to maintain SNR (almost) constant
        for i = 1:nW
            indx1 = (i-1)*wL + 1;
            indx2 = indx1 + wL - 1;
            pNoiseColored = sum(nColored(indx1:indx2).^2);
            pSignal = sum(y(indx1:indx2,1).^2);
            nColoredScaled  = nColored(indx1:indx2) * sqrt((pSignal/(10^(SNR/10)))/pNoiseColored);
            noise(indx1:indx2,1) = nColoredScaled;
        end 
    case 'experimental'
        load('noise_voluntary.mat')
        %-- Experimental noise has 213 realizations, each 60s long. To make it 120s, we repeat the same sequence.
        rNo1 = randi([1,213],1); %-- random experimental noise realization
        rNo2 = randi([1,213],1); %-- random experimental noise realization
        
        n1 = transpose(noise(rNo1,:));
        n2 = transpose(noise(rNo2,:));
        %-- Experimental noise is recorded at 1kHz but the required noise is at Fs (Hz) 
        noiseFs = 1/noiseTs;
        decimation = noiseFs/Fs;
        n1_d = decimate(nldat(n1,'domainIncr',noiseTs),decimation); n1_d = n1_d.data;
        n2_d = decimate(nldat(n2,'domainIncr',noiseTs),decimation); n2_d = n2_d.data;
        
        noiseNotScaled = [n1_d;n2_d];
        
        %-- Scaling the colored noise to maintain SNR (almost) constant
        noise = zeros(N,1);
        for i = 1:nW
            indx1 = (i-1)*wL + 1;
            indx2 = indx1 + wL - 1;
            pNoise = sum(noiseNotScaled(indx1:indx2).^2);
            pSignal = sum(y(indx1:indx2,1).^2);
            noiseScaled  = noiseNotScaled(indx1:indx2) * sqrt((pSignal/(10^(SNR/10)))/pNoise);
            noise(indx1:indx2,1) = noiseScaled;
        end         
end


