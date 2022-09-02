% function pvnlnblDemo
clearvars
close all

%% Add path to the NLID toolbox on your local system, where you have cloned or downloaded the NLID toolbox
% run 'S:\Biomed\REKLAB\myStuff\esobha1\RA 2021\Source Code\NPNPVH Latest NLID\initPath'

%% Load the identified model from Ehsan's implementation of the NPN-Hammerstein code
load('.\sim_data\example_PVHSimData_IEEEAccess2022.mat','udotD','schedVar','tq','Ts','decimation',...
                                                                                       'irf_len_r','nside_r','alfa',...
                                                                                       'q','m','n','p',...
                                                                                       'max_iter','threshold');
                                                                                   
%% First, create nldat data objects
velD = nldat(udotD,'domainIncr',Ts);
TQ = nldat(tq,'domainIncr',Ts); 
rho = nldat(schedVar,'domainIncr',Ts);

Fs = 1/Ts;

%% Measurement Noise
SNR = 15; %10; %15; %40 %100 
noiseType = 'tv-colored'; %'experimental'; %'tv-white'; 'tv-colored'; 'tiv'; 
wLength = 0.5;   %-- in seconds

noiseSource = 'new';

switch noiseSource 
    case 'new'
        switch noiseType 
            case 'tiv'
                noise = snr_gen(TQ.dataSet,SNR);
            case 'tv-white'
                noise = snr_tv_gen(TQ.dataSet,SNR,wLength/Ts,'white',Fs);
            case 'tv-colored'
                noise = snr_tv_gen(TQ.dataSet,SNR,wLength/Ts,'colored',Fs);
            case 'experimental'
                wLength = 1;
                noise = snr_tv_gen(TQ.dataSet,SNR,wLength/Ts,'experimental',Fs); 
        end
        noise = nldat(noise,'domainIncr',Ts);
    otherwise
        disp('There is no other noise source implemented')
end

%% Add noise to the output
TQ_n = TQ + noise;

%% Input/SV/Output (I/SV/O) Data structure 
z = cat(2,velD,TQ_n);  
set(z,'chanNames',{'perturbation velocity','reflex torque'},...
      'chanUnits',{'(rad/s)','Nm'});

sv = rho;
set(sv,'chanNames',{'SV, position'},...
      'chanUnits',{'(rad)'});
%% Instantiate a PV Hammerstein object from the pvnlbl class
pvhSys = pvnlbl;
ioNames = z.chanNames;
svNames = sv.chanNames;
set(pvhSys,'inputName',ioNames{1,1},'schedVarName',svNames{1,1},'outputName',ioNames{1,2});

%% Set identification method and its parameters
set(pvhSys,'idMethod','nppv-h');
%++ Set parameters
pvhSys.n = 7; %-- Default is 9
pvhSys.p = 3; %-- Default is 7

pvhSys.q = 8; %-- Default is 8
pvhSys.m = 3; %-- Default is 8

%% Identify the system using I/SV/O data structure
pvhSys = nlident(pvhSys,z,sv,'idMethod',pvhSys.idMethod,'decimation',decimation);

%% Simulate the identified PV-Hammerstein model
u = z(:,1);
u_d = decimate(u,decimation);
sv_d = decimate(sv,decimation);
TQ_d_hat = nlsim(pvhSys,u_d,sv_d);

%% Plot predicted output against measured output
TQ_d = decimate(TQ,decimation);
figure;
time = (0:length(TQ_d.dataSet)-1)*Ts*decimation;
plot(time,TQ_d.dataSet); hold on; plot(time,TQ_d_hat.dataSet,'r'); legend('Measured','Predicted')
xlabel('time (s)'); ylabel('torque (Nm)')
v = vaf(TQ_d,TQ_d_hat);
title(sprintf('Identification VAF = %0.1f%%',v.dataSet))

%% Plot the identified PV Hammerstein system
figure;
plot(pvhSys,'n_bins_input',50,'n_bins_sv',50); colormap('jet');

%% Plot static NL mimo basis
%++ Extracting SV Tchebychev polynomials from model.static_nl
PVNL = pvhSys.elements{1,1};           
figure;
plot(PVNL,'n_bins_input',50,'n_bins_sv',50); colormap('jet');

%% Plot dynamic IRF mimo basis
PVIRF = pvhSys.elements{1,2};  
figure;
plot(PVIRF,'n_bins_input',80,'n_bins_sv',50); colormap('jet');






