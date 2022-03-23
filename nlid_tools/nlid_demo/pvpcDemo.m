% function pvpcDemo
clearvars
close all

%% Add path to the NLID toolbox on your local system, where you have cloned or downloaded the NLID toolbox
run 'S:\Biomed\REKLAB\myStuff\esobha1\RA 2021\Source Code\NPNPVH Latest NLID\initPath'

%% Load data from a {PT,UT} trial pair of the pilot experimental data used for IEEE TBME publication
load('.\sim_data\IEEETBME2022_PilotData_AV_PT1_UT2.mat');  %% This contains input/output data in z (nldat object) and scheduling variable in rho (nldat object)

Ts = io.domainIncr;
Fs = 1/Ts;

input = io(:,1);
output = io(:,2);

%% Input/SV/Output (I/SV/O) Data structure 
z = cat(3,input,sv,output);
set(z,'chanNames',{'position perturbation','SV, ankle position','torque in response to perturbation'},...
      'chanUnits',{'(rad)','(rad)','Nm'});

%% Instantiate a PV-PC object from the pvpc class
pvpcStiffness = pvpc;
chanNames = z.chanNames;
set(pvpcStiffness,'inputName',chanNames{1,1},'schedVarName',chanNames{1,2},'outputName',chanNames{1,3});

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
pvpcStiffness = nlident(pvpcStiffness,z,'idMethod',pvpcStiffness.idMethod,'decimation',decimation,'rDelay',rDelay);

%% Plot output prediction against output
TQ_d = decimate_kian(output,decimation);
figure;
time = (0:length(TQ_d.dataSet)-1)*Ts*decimation;
subplot(4,1,1)
plot(time,TQ_d.dataSet); hold on; plot(time,pvpcStiffness.identTQt.dataSet,'r'); legend('Measured','Predicted')
title(sprintf('Identification VAF = %0.1f%%',pvpcStiffness.identVAF))
ylabel('Total torque (Nm)')
subplot(3,1,2)
plot(time,pvpcStiffness.identTQi.dataSet,'m')
ylabel('Intrinsic torque (Nm)')
subplot(3,1,3)
plot(time,pvpcStiffness.identTQr.dataSet,'k')
ylabel('Reflex torque (Nm)')
xlabel('Time (s)'); 

%% Plot the identified PV IRF model of intrinsic pathway ==> Still non-functional - TBD
% figure;
% PVIRF_i = pvpcStiffness.elements{1,1};
% plot(PVIRF_i,'n_bins_input',80,'n_bins_sv',50)

%% Plot the identified PV Hammerstein system of reflex pathway
figure;
PVH_r = pvpcStiffness.elements{2,1};
plot(PVH_r,'n_bins_input',50,'n_bins_sv',50)

