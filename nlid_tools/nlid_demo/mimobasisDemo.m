% function mimobasisDemo
clear all
close all
clear classes
clc

%% Add path to the NLID toolbox
addpath 'S:\Biomed\REKLAB\myStuff\esobha1\RA 2021\reklab_public\nlid_tools'
addpath 'S:\Biomed\REKLAB\myStuff\esobha1\RA 2021\reklab_public\utility_tools'
addpath 'S:\Biomed\REKLAB\myStuff\esobha1\RA 2021\reklab_public\stiffnessID'
addpath 'S:\Biomed\REKLAB\myStuff\esobha1\RA 2021\reklab_public\smi-2.0-devel'

%% Load the identified model from Ehsan's implementation of the NPN-Hammerstein code
load('S:\Biomed\REKLAB\myStuff\esobha1\RA 2021\Data Files\ehsanNPVHammCodeResults.mat','model','min_sv','max_sv','min_u','max_u');

%% First, test a static NL mimo basis
%++ Extracting SV Tchebychev polynomials from model.static_nl
sv_polynoms = cell(model.static_nl.inputExpOrder+1,1);

for j = 1:length(sv_polynoms)
    sv_polynoms{j,1} = polynom('polyCoef',model.static_nl.coeffs(j,:)',...
                               'polyOrder',model.static_nl.svExpOrder,...
                               'polyType','tcheb');                             
end

aa = mimobasis;
aa = set(aa,'coeffs',sv_polynoms,...
            'svRange',[min_sv,max_sv],...
            'inputRange',[min_u,max_u],...
            'svExpType','tcheb',...
            'inputExpType','tcheb',...
            'svExpOrder',model.static_nl.svExpOrder,...
            'inputExpOrder',model.static_nl.inputExpOrder);
           
figure;
plot(aa,'n_bins_input',50,'n_bins_sv',50)


%% Second, test a dynamic IRF mimo basis
%++ Extracting SV Tchebychev polynomials from model.static_nl
sv_polynoms = cell(model.h_r.LaguerreExpOrder+1,1);

for j = 1:length(sv_polynoms)
    sv_polynoms{j,1} = polynom('polyCoef',model.h_r.coeffs(j,:)',...
                               'polyOrder',model.h_r.svExpOrder,...
                               'polyType','tcheb');                             
end

min_u = 0;
max_u = (model.h_r.nLags - 1)*model.h_r.Ts;

bb = mimobasis;
bb = set(bb,'coeffs',sv_polynoms,...
            'svRange',[min_sv,max_sv],...
            'inputRange',[min_u,max_u],...
            'svExpType','tcheb',...
            'inputExpType','laguerre',...
            'svExpOrder',model.h_r.svExpOrder,...
            'inputExpOrder',model.h_r.LaguerreExpOrder,...
            'alfa',model.h_r.LaguerreAlfa);

figure;
plot(bb,'n_bins_input',80,'n_bins_sv',50)






