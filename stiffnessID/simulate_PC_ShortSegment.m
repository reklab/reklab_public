function [trqIntrinsicPredict,trqReflexPredict,trqTotalPredict,idVAF,posMeasured,trqMeasured] = simulate_PC_ShortSegment(intrinsic,reflex,z)
%This function simulates the PC system to segdat inputs
%It identifies initial conditions as part of the simulation
decimationRatio = 1;
zn=nldat(z)
position = zn(:,1);
position = position - mean(position,'omitnan');
torque = zn(:,2);
torque = torque - mean(torque,'omitnan');
torque = get(torque,'dataSet');
samplingTimeIntrinsic = get(intrinsic,'domainIncr');
samplingTimeReflex = get(reflex{2},'domainIncr');
if ~(samplingTimeReflex == samplingTimeIntrinsic)
    error('Intrinsic and reflex system sampling time mismatch')
else
    samplingTimeSystem = samplingTimeIntrinsic;
end
samplingTimeData = get(z,'domainIncr');
if (samplingTimeSystem < samplingTimeData)
    error('System sampling time must be equal or larger than the data sampling time')
elseif ~(samplingTimeData == samplingTimeSystem)
    decimationRatio = samplingTimeSystem / samplingTimeData;
    disp('System and data sampling time mismatch')
    disp(['Decimating Data by factor of : ',num2str(decimationRatio)]);
    samplingTime = samplingTimeSystem;
end

segmentOnsetPointer = get(z,'onsetPointer');
inputSegmentOnsetPointer = segmentOnsetPointer (:);
outputSegmentOnsetPointer = segmentOnsetPointer (:);
segmentLength = get(z,'segLength');
inputSegmentLength = segmentLength (:);
outputSegmentLength = segmentLength (:);
if ~( isequal(inputSegmentOnsetPointer,outputSegmentOnsetPointer) &&...
        isequal(inputSegmentLength,outputSegmentLength))
    error('The input and output segment onset pointers and lengths must be equal')
else
    segmentOnsetPointer = segmentOnsetPointer(:);
    segmentLength = segmentLength(:);
	segmentEndpointer = segmentOnsetPointer + segmentLength - 1;

end
intrinsicIRFCoeff = samplingTime * get(intrinsic,'dataSet');

intrinsicIRF_Length = abs(intrinsic.domainStart / samplingTime);
lagsIntrinsic = (-intrinsicIRF_Length:1:intrinsicIRF_Length);
numLagsIntrinsic = length(intrinsic.dataSet);
positionDelay = zeros(size(position,1),numLagsIntrinsic);
for j = 1 : numLagsIntrinsic
	posDelay = del(position,lagsIntrinsic(j) * samplingTime);
	positionDelay(:,j) = get(posDelay,'dataSet');
end
positionNldat = nldat(get(position,'dataSet'),'domainIncr',samplingTimeData);
velocity = ddt(positionNldat);
velocity = get(velocity,'dataSet');
    
pointer = 1;
switch_time = zeros(length(segmentEndpointer)-1,1);
positionDelaySegments = zeros(sum(segmentLength),numLagsIntrinsic);
velocitySegments = zeros(sum(segmentLength),1);
trqMeasured = zeros(sum(segmentLength),1);
for i = 1 : length(segmentEndpointer)
    vel_seg = velocity(segmentOnsetPointer(i):segmentEndpointer(i));
    velocitySegments(pointer:pointer+segmentLength(i)-1) = vel_seg;
    trqMeasured(pointer:pointer+segmentLength(i)-1) = torque(segmentOnsetPointer(i):segmentEndpointer(i));
    positionDelaySegments(pointer:pointer+segmentLength(i)-1,:) = positionDelay(segmentOnsetPointer(i):segmentEndpointer(i),:);
    pointer = pointer + segmentLength(i);
    switch_time(i) = pointer;
end
[velocitySegments,~,~,~] = decimate_segment(velocitySegments,switch_time(1:end-1),decimationRatio);
positionDelaySegments = bsxfun(@minus,positionDelaySegments,mean(positionDelaySegments));
u_i = zeros(size(velocitySegments,1),numLagsIntrinsic);
for i = 1 : numLagsIntrinsic
	[u_i(:,i),~,~,~] = decimate_segment(positionDelaySegments(:,i),switch_time(1:end-1),decimationRatio);
end
posMeasured = u_i(:,ceil(numLagsIntrinsic/2));
[trqMeasured,switch_time,segLength,~] = decimate_segment(trqMeasured,switch_time(1:end-1),decimationRatio);
velocitySegments = segdat(velocitySegments,'onsetPointer',switch_time(1:end-1),'segLength',segLength,'domainIncr',samplingTime);
posMeasured = segdat(posMeasured,'onsetPointer',switch_time(1:end-1),'segLength',segLength,'domainIncr',samplingTime);
trqMeasured = segdat(trqMeasured,'onsetPointer',switch_time(1:end-1),'segLength',segLength,'domainIncr',samplingTime);
trqIntrinsicPredict = u_i * intrinsicIRFCoeff;
tqIResidual = trqMeasured - trqIntrinsicPredict;
trqIntrinsicPredict = segdat(trqIntrinsicPredict,'domainIncr',samplingTime,'chanNames','Torque (Nm)','comment','Intrinsic Torque','onsetPointer',switch_time(1:end-1),'segLength',segLength);
tqIResidual = segdat(tqIResidual,'domainIncr',samplingTime,'chanNames','Torque (Nm)','comment','Intrinsic Torque','onsetPointer',switch_time(1:end-1),'segLength',segLength);
zReflex = cat(2,velocitySegments,tqIResidual);
zReflex.domainStart=(zReflex.onsetPointer-1)*zReflex.domainIncr;
trqReflexPredict = nlsim(reflex,zReflex);
trqIntrinsicPredict.domainStart=trqReflexPredict.domainStart;
trqMeasured.domainStart=trqReflexPredict.domainStart;
trqTotalPredict = trqIntrinsicPredict + trqReflexPredict;
idVAF = vaf(trqMeasured,trqTotalPredict);
end
function [output,switch_time_new,interval,p] = decimate_segment(input,switch_time,decimation_ratio)
    sw = [1;switch_time];
    p = length(sw);
    sw = [sw; length(input)+1];
    output = [];
    switch_time_new = 1;
    interval=zeros(length(switch_time),1);
    segment_onset=1;
    for i = 1 : p
        output_temp = decimate(input(sw(i):sw(i+1)-1),decimation_ratio);
        output = [output;output_temp];
        segment_onset = segment_onset+length(output_temp);
        switch_time_new = [switch_time_new;segment_onset];
        interval(i) = length(output_temp);
    end
end