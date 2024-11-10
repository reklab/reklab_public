%% This is a code to that takes the Hammerstein model identified between 
%% Uterine Activity (UA) and Fetal Heart Rate (FHR) using NLID toolbox (chebychev NL + state-space model), 
%% and converts the chebyhev static NL to a half-wave-rectifier (HWR) nonlinearity.
%% The code reads the identified from a MAT file that REK shared with me on March 9, 2023, which has 15 trials,
%% but only the following trials are unique: [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 15].  

clear
close all
clc

%%
set(0,'DefaultFigureWindowStyle','docked');

uniqueTrials = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 15];
results = cell(length(uniqueTrials),1);

Ts = 2;

Ti = 0;
Tf = 300;
Tvec = Ti:Ts:Tf;

nSlopes = 100;      %50; 100; 150; 200;
nThresholds = 100;  %50; 100; 150; 200;

%% Load the identified models as a table
matFilesPath = [pwd,'\Data\'];

%-- NLID table
dataFile = 'snl2hwrDemo.mat';  %-- Data are in a table named TV
load(dataFile);
TV_NLID = TV;

%% Comparison of the predictions
for i = 1:length(uniqueTrials)
    trialNo = uniqueTrials(i);
    %-- Read input
    u = TV_NLID.data{trialNo,1}(:,1);
    y = TV_NLID.data{trialNo,1}(:,2);
    
    % Delay the input
    ud = del(u,TV_NLID.delay(trialNo,1)*Ts);
    
    %====================================
    %% Comparing the identified systems
    %====================================
    %% Normalize NLID Models
    M_NLID = TV_NLID.model{trialNo,1};
    NL_NLID = M_NLID.elements{1,1};
    LD_NLID = M_NLID.elements{1,2};
    
    %-- Normalize the gain of the linear dynamics to 1
    LD_NLID_SS = ss(LD_NLID);
    dcGain_NLID = dcgain(LD_NLID_SS);
    LD_NLID_SS_N = LD_NLID_SS / dcGain_NLID;

    LD_NLID_N = ssm('A',LD_NLID_SS_N.A,'B',LD_NLID_SS_N.B,'C',LD_NLID_SS_N.C,'D',LD_NLID_SS_N.D,'domainIncr',LD_NLID_SS_N.Ts,'orderSelect','manual');

    %-- Apply the normalization to static NL
    NL_NLID_N = NL_NLID;
    NL_NLID_N.polyCoef = NL_NLID.polyCoef * dcGain_NLID;

    M_NLID_N = M_NLID;
    M_NLID_N.elements = {NL_NLID_N, LD_NLID_N};
    
    %% Finding HWR fits using curve fit optimization criterion
    tic; [HWR_c,info_c] = snl2hwr(NL_NLID_N,ud,'optim_criterion','curve_err','slope_range',[-10,10],'plot_figures',true,'n_slopes',nSlopes,'n_thresholds',nThresholds); toc
    
    %% Finding HWR fits using prediction fit to NL output optimization criterion
    tic; [HWR_p,info_p] = snl2hwr(NL_NLID_N,ud,'optim_criterion','pred_err','slope_range',[-10,10],'plot_figures',true,'n_slopes',nSlopes,'n_thresholds',nThresholds); toc

end
