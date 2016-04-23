clc
%loading experimental data
load experimental_data.mat
data = z_pf;
%select the position and torque from z_pf. You can try other data records
%Each variable has position and torque records as input and output signals
position = get(data(:,1),'dataSet');%input position
torque = get(data(:,2),'dataSet');%output torque
samplingTime = get(data,'domainIncr');
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
disp('Identification finished')
disp(['Identification VAF was : ',num2str(vafs(1))])
%% Plotting System
figure
subplot(2,1,1)
zIntrinsic = cat(2,decimate(data(:,1),10),nlsim(intrinsic, decimate(data(:,1),10)));
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
pause
%% Validation
% The purpose of this section is to simulate the response of 
% the identified system to novel segments
% Generating new segments from the data
disp('Validation Starting')
numSegmentValidation = 20;
minSegment = 0.25;
maxSegment = 0.5;
minSegment = floor(minSegment / samplingTime);
maxSegment = floor(maxSegment / samplingTime);

disp(['Randomly selecting : ',num2str(numSegmentValidation),' segments for validation'])
segLengthValidation = randi([minSegment,maxSegment],numSegmentValidation,1);% vector of segment lengths
onsetPointerValidation = randi([1,length(position) - maxSegment],numSegmentValidation,1);%vector of segment onset
positionValidation = segdat(position,'onsetPointer',onsetPointerValidation,'segLength',segLengthValidation,'domainIncr',samplingTime...
,'comment','Position','chanNames','Joint angular position (rad)');
torqueValidation = segdat(torque,'onsetPointer',onsetPointerValidation,'segLength',segLengthValidation,'domainIncr',0.001...
,'comment','Torque','chanNames','Joint torque (Nm)');
%Concatenate input and output
zValidation = cat(2,positionValidation,torqueValidation);
[tqIValidation,tqRValidation,tqTValidation,vafValidation,posMeasured,trqMeasured] = simulate_PC_ShortSegment(intrinsic,reflex,zValidation);
disp(['Validation VAF was : ',num2str(vafValidation(1)),' for validation'])
%Let's visualize
onsetPointer = tqIValidation.onsetPointer;
segLength = tqIValidation.segLength;
posMeasured = nldat(posMeasured.dataSet,'domainIncr',0.01,'comment',['Position; VAF = ',num2str(vafValidation(1))],'chanNames','Position (rad)');
trqMeasured = nldat(trqMeasured.dataSet,'domainIncr',0.01,'comment','Total Torque','chanNames','Torque (Nm)');
tqIValidation = nldat(tqIValidation.dataSet,'domainIncr',0.01,'comment','Intrinsic Torque','chanNames','Torque (Nm)');
tqRValidation = nldat(tqRValidation.dataSet,'domainIncr',0.01,'comment','Reflex Torque','chanNames','Torque (Nm)');
tqTValidation = nldat(tqTValidation.dataSet,'domainIncr',0.01,'comment','Total Torque','chanNames','Torque (Nm)');

disp('Plotting validation segments')
for i = 1 : numSegmentValidation
	h = figure;
    subplot(4,1,1)
    plot(posMeasured(onsetPointer(i):onsetPointer(i)+segLength(i)-1,1))
    subplot(4,1,2)
    plot(trqMeasured(onsetPointer(i):onsetPointer(i)+segLength(i)-1,1))
    hold on
    plot(tqTValidation(onsetPointer(i):onsetPointer(i)+segLength(i)-1,1),'line_color','r')
    subplot(4,1,3)
    plot(tqIValidation(onsetPointer(i):onsetPointer(i)+segLength(i)-1,1),'line_color','r')
    subplot(4,1,4)
    plot(tqRValidation(onsetPointer(i):onsetPointer(i)+segLength(i)-1,1),'line_color','r')
    pause
    close(h)
end