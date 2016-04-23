clc
%loading experimental data
load experimental_data.mat
%select the position and torque from z_pf. You can try other data records
%Each variable has position and torque records as input and output signals
position = get(z_pf(:,1),'dataSet');%input position
torque = get(z_pf(:,2),'dataSet');%output torque
samplingTime = get(z_pf,'domainIncr');
%Randomly select segments from this data
minSegment = 0.5;%minimum segment length in s
maxSegment = 1;%maximum segment length in s
numSegment = 50;%number of segments
minSegment = floor(minSegment / samplingTime);
maxSegment = floor(maxSegment / samplingTime);
segLength = randi([minSegment,maxSegment],numSegment,1);% vector of segment lengths
onsetPointer = randi([1,length(position) - maxSegment],numSegment,1);%vector of segment onset
%Define segdat objects
%segdat is a child of nldat (data class of NLID toolbox) for segmented data
position = segdat(position,'onsetPointer',onsetPointer,'segLength',segLength,'domainIncr',samplingTime...
,'comment','Position','chanNames','Joint angular position (rad)');
torque = segdat(torque,'onsetPointer',onsetPointer,'segLength',segLength,'domainIncr',0.001...
,'comment','Torque','chanNames','Joint torque (Nm)');
%Concatenate input and output
z = cat(2,position,torque);
%Let's visualize the data
h = figure;
subplot(2,1,1)
plot(position)
subplot(2,1,2)
plot(torque)
xAxisPanZoom
disp('Press any key to continue')
pause
if (ishandle(h))
    close(h);
end
%% Identification of the PC structure
% Identify the system. intrinsic and reflex are nlid objects
[intrinsic, reflex, tqI, tqR, tqT, vafs] = SS_SDSS_stiffnessID (z);
%% Plotting System
figure
subplot(2,1,1)
zIntrinsic = cat(2,decimate(z_pf(:,1),10),nlsim(intrinsic, decimate(z_pf(:,1),10)));
frespIntrinsic = fresp(zIntrinsic);
FrespData = frespIntrinsic.dataSet;
F = nldat(abs(FrespData(:,1)),'domainIncr',get(frespIntrinsic,'domainIncr'));
plot(F)
set(gca,'XScale','log');
xlim([0.1,50])
title(['Intrinsic FRF; VAF = ',num2str(vafs(1))])
ylabel('Magnitude (dB)')
xlabel('Frequency (Hz)')
subplot(2,2,3)
plot(reflex{1})
title('Reflex Static Nonlinearity')
xlim([-2,2])
subplot(2,2,4)
F = fresp(reflex{2});
F = nldat(abs(F.dataSet),'domainIncr',get(F,'domainIncr'));
plot(F)
set(gca,'XScale','log');
title('Reflex Linear Dynamics FRF')
xlim([0.1,50])
ylabel('Magnitude (dB)')
xlabel('Frequency (Hz)')
%% Validation
% The purpose of this section is to simulate the response of 
% the identified system to novel segments