% function pvnlDemo
clearvars
close all

%% Add path to the NLID toolbox
run 'S:\Biomed\REKLAB\myStuff\esobha1\RA 2021\Source Code\NPNPVH Latest NLID\initPath'

%% Load the identified model from Ehsan's implementation of the NPN-Hammerstein code
load('S:\Biomed\REKLAB\myStuff\esobha1\RA 2021\Data Files\ehsanNPVHammCodeResults.mat','udotD','schedVar','tq','TQhat','Ts','decimation','model');

%% First, create nldat data objects
u = nldat(udotD,'domainIncr',Ts);
y = nldat(tq,'domainIncr',Ts); 
rho = nldat(schedVar,'domainIncr',Ts);

%% Decimate the signals
u_d = decimate_kian(u,decimation);
y_d = decimate_kian(y,decimation);
rho_d = decimate_kian(rho,decimation);

%% Simulate the identified NPNPV-H model captured in the the "model" structure
[yp,zp1] = pvHammLaguerreSim(model,u_d,rho_d);

delta = yp.dataSet - TQhat.Data;
v = vaf(yp.dataSet,TQhat.Data);
figure; 
subplot(2,1,1); plot(delta); 
subplot(2,1,2); plot(yp.dataSet,'r'); hold on; plot(TQhat.Data,'b'); legend('Simulated','Identified') 

%% Simulate the identified static PVNL element of the identified NPNPV-H model
zp2 = pvnlSim(model.static_nl,u_d,rho_d);

deltaZ = zp1 - zp2;
figure; 
subplot(2,1,1); plot(deltaZ); 
subplot(2,1,2); plot(zp1,'lineColor','r'); hold on; plot(zp2,'lineColor','b'); legend('pvHammLaguerreSim','pvnlSim')
