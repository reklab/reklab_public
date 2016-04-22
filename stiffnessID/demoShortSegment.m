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
numSegment = 60;%number of segments
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
%%
system = SS_SDSS_stiffnessID (z);
figure
subplot(2,1,1)
plot(system{1})
title('Intrinsic IRF')
subplot(2,2,3)
plot(system{2}{1})
subplot(2,2,4)
plot(system{2}{2})
