clear all
close all
clc

%% Specifying random duration specifications for the Position Perturbation Signal Pulse Widths
minT = 160; %-- in milliseconds or in samples (i.e., with 1msec sampling = 160 msec)
maxT = 190; %-- in milliseconds or in samples (i.e., with 1msec sampling = 190 msec)

%% Number od perturbation realizations
Nr = 1; %10; %2 %200; %40;  % 33%--- Number of realizations

%% Sampling time and record length
Ts = 0.001;  %-- sampling time in sec
Nxdot = 120000; %120000; %--- Total number of time samples in the perturbation sequence  (120000 samples = 120 seconds was used in data Generated for Chris on Aug. 10, 2012)

%% Time signal for plotting
time = (0:Ts:(Nxdot-1)*Ts)';

%% Switching and Scale Factor of Velocity
nSwitch = 1; %1; %2; %5; %4; %3;   %-- Switching rate between Pos and Neg velocities (in samples)
desiredpdfType = 'exponential_2';  %'V-shape','exponential_1','exponential_2';

Npdf = 100000;          %--- Total number of samples in the realization of the desired PDF
Nbins = 200;

[xPertData,xdotPertData] = gen_pralds(minT,maxT,nSwitch,Ts,Nxdot,Nr,desiredpdfType,Npdf,Nbins,'');

figure(102); 
subplot(2,1,1); 
plot(time,xPertData,'b')
xlabel('Time (s)')
ylabel('PRALDS Amplitude (rad)')
subplot(2,1,2)
plot(time,xdotPertData,'b')
xlabel('Time (s)')
ylabel('PRALDS Velocity (rad/s)')

figure;
hist(xdotPertData,Nbins)
xlabel('PRALDS Velocity (rad/s)')
ylim([0,100])
ylabel('Number of occurrence')