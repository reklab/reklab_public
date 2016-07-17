%%% Demo for the Structural Decomposition and SubSpace (SDSS) method
clear
clc
close all
disp('Identifiying PC structure; experimental data PF condition')
load experimental_data.mat
z=z_pf2;
[intrinsic, reflex, tqI, tqR, tqT] = SDSS_stiffnessID (z);
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
V = vaf(decimate(z(:,2),10),tqT);
disp(['VAF was: ',num2str(V.dataSet)])
%%
disp('Identifiying PC structure; experimental data DF condition')
z=z_df3;
[intrinsic, reflex, tqI, tqR, tqT] = SDSS_stiffnessID (z,'delay',0.05);
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
V = vaf(decimate(z(:,2),10),tqT);
disp(['VAF was: ',num2str(V.dataSet)])