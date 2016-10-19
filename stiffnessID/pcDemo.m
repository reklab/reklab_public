%% pcDemo
load experimental_data.mat
z=z_pf2;
[hIntrinsic, mReflex, nlmStiff, VAFT, VAFi, VAFr, tqI, tqR, tqT]= PC_stiffnessID (z);
%% Yong's state-space method
[system_ss, pIntrinsic, tqPredicted]= SS_stiffnessID (z);
%% SDSS state-space method
[intrinsic, reflex, tqI, tqR, tqT] = SDSS_stiffnessID (z);
%%
figNum=1;
xLim=[10 20];
pos=decimate(z(:,1),1);
tq=decimate(z(:,2),10);
vafI=vaf(tq,tqI);
vafR=vaf(tq-tqI,tqR);
vafT=vaf(tq,tqT)

figure(figNum);
subplot (4,1,1);
plot(pos);set(gca,'xlim',xLim),set(gca,'xlim',xLim,'xticklabel',{}); 
; 
title('Position');
xlabel('');
ylabel('rad');


subplot (4,1,2);
plot(tq,'line_color','r');
set(gca,'xlim',xLim,'xticklabel',{}); 
title('Torque');
xlabel('');
ylabel('Nm');
set(gca,'xlim',xLim,'xticklabel',{}); 

subplot (4,1,3);
title('GS EMG');
xlabel('');
set(gca,'xlim',xLim,'xticklabel',{}); 


set(gca,'xlim',xLim); 
subplot (4,1,4);
title('TA EMG')
xlabel('');

set(gca,'xlim',xLim); 
fig_mod(figNum);

figNum=figNum+1;
figure(figNum);  clf
subplot (4,1,1);
plot(pos); set(gca,'xlim',xLim, 'xticklabel',{}); 
title('Position');
xlabel('');
ylabel('rad');
subplot (4,1,2); plot(tqI,'line_color','m');
set(gca,'xlim',xLim,'xticklabel',{}); 
title('Intrinsic Torque')
ylabel('Nm')
xlabel('');
subplot  (4,1,3); plot (tqR, 'line_color','g');
set(gca,'xlim',xLim, 'xticklabel',{}); 
title('Reflex Torque')
ylabel('Nm');
xlabel('');
subplot (4,1,4);
plot(tq,'line_color','r');
h=line(tqT+mean(tq));
set(h,'color','cyan')
set(gca,'xlim',xLim); 
title('Total Torque');
ylabel('Nm');
legend('Observed','Predicted','Location','Southwest');
fig_mod(figNum);




