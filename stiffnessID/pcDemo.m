function pcDemo
%% Demonstrate various aproaches to the identificatio of the parallel cascad emode 
load experimental_data.mat
z=z_pf2;
%% Rectifier followed by IRF 
reflexIdMethod='rect'; 
[hIntrinsic, mReflex, nlmStiff, VAFT, VAFi, VAFr, tqI, tqR, tqT]= PC_stiffnessID (z,'reflex_id_method',reflexIdMethod);
figNum=1;
pcDemoPlot (z, tqI, tqR, tqT, nlmStiff, figNum)


%%  Estimate reflex pathway using HK method 
reflexIdMethod='hk'; 
[hIntrinsic, mReflex, nlmStiff, VAFT, VAFi, VAFr, tqI, tqR, tqT]= PC_stiffnessID (z,'reflex_id_method',reflexIdMethod);
figNum=3;
pcDemoPlot (z, tqI, tqR, tqT, nlmStiff, figNum)

%%  Estimate reflex pathway using sls method 
reflexIdMethod=''; 
[hIntrinsic, mReflex, nlmStiff, VAFT, VAFi, VAFr, tqI, tqR, tqT]= PC_stiffnessID (z,'reflex_id_method',reflexIdMethod);
figNum=5;
pcDemoPlot (z, tqI, tqR, tqT, nlmStiff, figNum)


%% Yong's state-space method
%%[system_ss, pIntrinsic, tqPredicted]= SS_stiffnessID (z);


%% SDSS state-space method
figNum=7;

[intrinsic, reflex, tqI, tqR, tqT] = SDSS_stiffnessID (z);
nlmStiff=nlm;
nlmStiff{1,1}=intrinsic;
nlmStiff{2,1}=reflex; 
figNum=5;
pcDemoPlot (z, tqI, tqR, tqT, nlmStiff, figNum)



%%
end
function pcDemoPlot (z, tqI, tqR, tqT, nlmStiff, figNum)
figNum=1;
xLim=[10 20];
pos=decimate(z(:,1),1);
tq=decimate(z(:,2),10);
vafI=vaf(tq,tqI);
vafR=vaf(tq-tqI,tqR);
vafT=vaf(tq,tqT)

figure(figNum);
subplot (2,2,1);
plot (nlmStiff{1,1}); 
subplot (2,2,3);
R=(nlmStiff{2,1})
plot(R{1,1});
subplot (2,2,4);
plot (R{1,2})
streamer (nlmStiff.comment); 


figMod(figNum);
%%
figNum=figNum+1;
figure(figNum);  clf
subplot (4,1,1);
plot(pos); set(gca,'xlim',xLim, 'xticklabel',{});
title('Position');
xlabel('');
ylabel('rad');
subplot (4,1,2); plot(tqI,'linecolor','m');
set(gca,'xlim',xLim,'xticklabel',{});
title('Intrinsic Torque')
ylabel('Nm')
xlabel('');
subplot  (4,1,3); plot (tqR, 'linecolor','g');
set(gca,'xlim',xLim, 'xticklabel',{});
title('Reflex Torque')
ylabel('Nm');
xlabel('');
subplot (4,1,4);
plot(tq,'linecolor','r');
h=line(tqT+mean(tq));
set(h,'color','cyan')
set(gca,'xlim',xLim);
title(['Total Torque %VAF=' num2str(vafT.dataSet)]);
ylabel('Nm');
legend('Observed','Predicted','Location','Southwest');
streamer(nlmStiff.comment); 
figMod(figNum);
end





