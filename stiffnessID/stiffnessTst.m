%% Test stiffness identification methods
load experimental_data.mat
z=z_pf2;
[hIntrinsic, mReflex, nlmStiff, VAFT, VAFi, VAFr, tqI, tqR, tqT]= PC_stiffnessID (z);
%% Yong's state-space method
[system_ss, pIntrinsic, tqPredicted]= SS_stiffnessID (z);
%% SDSS state-space method
[intrinsic, reflex, tqI, tqR, tqT] = SDSS_stiffnessID (z);
