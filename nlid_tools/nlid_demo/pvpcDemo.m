% function pvpcDemo
clearvars
close all

%% Add path to the NLID toolbox on your local system, where you have cloned or downloaded the NLID toolbox
run 'S:\Biomed\REKLAB\myStuff\esobha1\RA 2021\Source Code\NPNPVH Latest NLID\initPath'

%% Load data from a {PT,UT} trial pair of the pilot experimental data used for IEEE TBME publication
load('.\sim_data\IEEETBME2022_PilotData_AV_PT1_UT2.mat','io','sv');  %% This contains input/output data in z (nldat object) and scheduling variable in rho (nldat object)

Ts = io.domainIncr;
Fs = 1/Ts;

input = io(:,1);
output = io(:,2);

%% Input/SV/Output (I/SV/O) Data structure 
z = io;  
set(z,'chanNames',{'position perturbation','torque in response to perturbation'},...
      'chanUnits',{'(rad)','Nm'});

set(sv,'chanNames',{'SV, ankle position'},...
      'chanUnits',{'(rad)'});

%% Instantiate a PV-PC object from the pvpc class
pvpcStiffness = pvpc;
ioNames = z.chanNames;
svNames = sv.chanNames;
set(pvpcStiffness,'inputName',ioNames{1,1},'schedVarName',svNames{1,1},'outputName',ioNames{1,2});

%% Set identification method and its parameters
set(pvpcStiffness,'idMethod','nppv-pc');
%++ Set the identification parameters
%== INTRINSIC Pathway
pvpcStiffness.irf_len_i = 0.04;  %-- in seconds
% pvpcStiffness.nside_i = 2;       %-- Number of sides of intrinsic IRF - Setting this gives an error! To be investigated. For now, default is 2.
pvpcStiffness.p_i = 2;           %0; 2; 3;

%== REFLEX Pathway
% pvpcStiffness.nside_r = 1;       %-- Number of sides of reflex IRF - Setting this gives an error! To be investigated. For now, default is 1.
rDelay = 0.05;         %-- in seconds
pvpcStiffness.irf_len_r = 0.8 - rDelay;
pvpcStiffness.alfa = 0.5;

pvpcStiffness.n = 7;   %-- Order of expansion of reflex NL w.r.t. delayed velocity. Default is 7.
pvpcStiffness.p = 2+1; %-- Order of expansion of reflex NL w.r.t. SV.               Default is 9.

pvpcStiffness.q = 8;   %-- Laguerre expansion order for reflex dynamics.    Default is 8.
pvpcStiffness.m = 0;   %-- Order of expansion of reflex dynamics w.r.t. SV. Default is 8.

pvpcStiffness.max_iter = 500; 
pvpcStiffness.threshold = 10^-10;                  %-- Threshold on SSE for terminitiaon of the iterative search'

%% Identify the system using I/SV/O data structure
decimation = 10;
pvpcStiffness = nlident(pvpcStiffness,z,sv,'idMethod',pvpcStiffness.idMethod,'decimation',decimation,'rDelay',rDelay);

%% Simulate the identified PV-Hammerstein model
u = z(:,1);
u_d = decimate_kian(u,decimation);
sv_d = decimate_kian(sv,decimation);
[tqT_d_hat, tqI_d, tqR_d] = nlsim(pvpcStiffness,u_d,sv_d,'rDelay',rDelay);

%% Plot predicted output against measured output, as well as the intrinsic and reflex torques
tqT_d = decimate_kian(output,decimation);
time = (0:length(tqT_d.dataSet)-1)*Ts*decimation;

figure;
subplot(5,1,1)
plot(u_d)
ylabel('input (rad)')
title('Position perturbation')

subplot(5,1,2)
plot(sv_d)
ylabel('SV (rad)')
title('Ankle position')

subplot(5,1,3)
plot(time,tqT_d.dataSet - mean(tqT_d.dataSet)); hold on; plot(time,tqT_d_hat.dataSet - mean(tqT_d_hat.dataSet),'r'); legend('Measured','Predicted')
ylabel('torque (Nm)')
v = vaf(tqT_d,tqT_d_hat);
title(sprintf('Total Torque Simulation VAF = %0.1f%%',v.dataSet))

subplot(5,1,4)
plot(tqI_d)
ylabel('torque (Nm)')
v = vaf(tqT_d,tqI_d);
title(sprintf('Intrinsic torque contribution to total torque, VAF = %0.1f%%',v.dataSet))

subplot(5,1,5)
plot(tqR_d)
ylabel('torque (Nm)')
v = vaf(tqT_d,tqR_d);
title(sprintf('Reflex torque contribution to total torque, VAF = %0.1f%%',v.dataSet))
xlabel('time (s)');

%% Plot the identified PV Hammerstein system of reflex pathway
figure;
PVH_r = pvpcStiffness.elements{2,1};
plot(PVH_r,'n_bins_input',50,'n_bins_sv',50)

