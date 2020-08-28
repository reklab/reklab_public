%%% Demo for the Structural Decomposition and SubSpace (SDSS) method
clear
clc
close all
disp('SDSS demo')
disp('Identifiying PC structure during tonic isometric contraction when the ankle plantarflexor muscles are active')
load experimental_data.mat
z=z_pf2;
[intrinsic, reflex, tqI, tqR, tqT] = SDSS_stiffnessID (z,'orderselectmethod',2);
figure
subplot(2,1,1)
plot(intrinsic)
title('IRF of intrinsic pathway')
subplot(2,2,3)
plot(reflex{1})
title('Reflex static nonlinearity')
subplot(2,2,4)
plot(reflex{2})
title('Reflex linear dynamics')
torqueDecimated = decimate(z(:,2),10);
positionDecimated = decimate(z(:,1),10);
torqueDecimated = torqueDecimated - mean(torqueDecimated);
set(torqueDecimated,'domainStart',0);
V = vaf(torqueDecimated,tqT);
disp(['VAF was: ',num2str(V.dataSet)])
figure
subplot(4,1,1)
plot(positionDecimated)
title('Position (rad)')
xlabel('')
ylabel('position (rad)')
subplot(4,1,2)
plot(tqI)
title('Intrinsic torque (Nm)')
xlabel('')
ylabel('Torque (Nm)')
subplot(4,1,3)
plot(tqR)
title('Reflex torque (Nm)')
xlabel('')
ylabel('Torque (Nm)')
subplot(4,1,4)
plot(torqueDecimated)
hold on
plot(tqT,'lineColor','r')
title('Total torque (Nm)')
xlabel('Time (s)')
ylabel('Torque (Nm)')
legend('Measured','Predicted')
xAxisPanZoom
disp('Use Zoom/Pan on the figure.')
disp('Press any key to continue.')
pause
%%
disp('Identifiying PC structure during tonic isometric contraction when the ankle dorsiflexor muscles are active')
z=z_df3;
[intrinsic, reflex, tqI, tqR, tqT] = SDSS_stiffnessID (z,'delay',0.05,'orderselectmethod',2);
figure
subplot(2,1,1)
plot(intrinsic)
title('IRF of intrinsic pathway')
subplot(2,2,3)
plot(reflex{1})
title('Reflex static nonlinearity')
subplot(2,2,4)
plot(reflex{2})
title('IRF of Reflex linear dynamics')
torqueDecimated = decimate(z(:,2),10);
positionDecimated = decimate(z(:,1),10);
torqueDecimated = torqueDecimated - mean(torqueDecimated);

V = vaf(decimate(z(:,2),10),tqT);
disp(['VAF was: ',num2str(V.dataSet)])
figure
subplot(4,1,1)
plot(positionDecimated)
title('Position (rad)')
xlabel('')
ylabel('position (rad)')
subplot(4,1,2)
plot(tqI)
title('Intrinsic torque (Nm)')
xlabel('')
ylabel('Torque (Nm)')
subplot(4,1,3)
plot(tqR)
title('Reflex torque (Nm)')
xlabel('')
ylabel('Torque (Nm)')
subplot(4,1,4)
plot(torqueDecimated)
hold on
plot(tqT,'lineColor','r')
title('Total torque (Nm)')
xlabel('Time (s)')
ylabel('Torque (Nm)')
legend('Measured','Predicted')
disp('Use Zoom/Pan on the figure.')
disp('Press any key to continue.')
xAxisPanZoom