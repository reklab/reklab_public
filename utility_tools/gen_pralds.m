function [xPertData,xdotPertData] = gen_pralds(minT,maxT,nSwitch,Ts,Nxdot,Nr,desiredpdfType,Npdf,Nbins,fileName)
%% This function generates a desired number of realizations from a position/velocity
%% perturbation signal that me and Kian have recently designed

%% Perturbation parameters
vMin = 0;
vMax = +1; %+1; %+2; %+3; %+1.5;
delta_v = 0.01;

switch desiredpdfType
    case 'V-shape'
        xdotScale = 12;  %8(pos: +/-0.015, vel: +/-1.26); %10(pos: +/-0.018, vel: +/-1.58); %12(pos: +/-0.023, vel: +/-1.9); %20; %-- Scale factor for velocity
    case 'exponential_2'
        xdotScale = 11;  %11(pos: +/-0.022, vel: +/-1.75);
end

%% Time signal for plotting
time = (0:Ts:(Nxdot-1)*Ts)';

%% Generating Nr realizations of the perturbation
xPertData = zeros(Nr,Nxdot);
xdotPertData = zeros(Nr,Nxdot);

%+++ We add the perturbations to an flb file cases, to Nr MAT files, and in a Data matrix that would be saved into a single MAT file 
h = waitbar(0,'Please wait for all the realizations of the PRALDS perturbation sequence to generate...');

for rr = 1:Nr
    waitbar(rr/Nr,h);
    
    %--- This is the code that generates one realization of the new perturbation
    [xdot,x] = newPert(vMin,vMax,delta_v,xdotScale,minT,maxT,nSwitch,desiredpdfType,Npdf,Nbins,Nxdot,Ts);
    %% writing Nr perturbations to Nr MAT file
    % fileName = ['xPert',num2str(rr),'.mat'];
    % eval(['save ',fileName,' x xdot vMin vMax minT maxT nSwitch desiredpdfType'])
    
    %% writing Nr perturbations into Nr cases in an FLB file
    % D.Data = x/max(abs(x));  %-- Data is normalized to 1,-1 
    % mat2flb('posPert_280212', 'append', D);
    
    %% writing Nr perturbations into a matrix and saving the results into a MAT file
    xPertData(rr,:) = x/max(abs(x)); 
    xdotPertData(rr,:) = xdot/max(abs(xdot)); 
end

% figure(101);
% plot(time,xPertData)
% xlabel('Time (s)')
% ylabel('PRALDS Amplitude')

close(h)
if ~isempty(fileName)
    save(fileName,'xPertData','xdotPertData','vMin','vMax','delta_v','xdotScale','minT','maxT','nSwitch','desiredpdfType','Npdf','Nbins','Nxdot')
end