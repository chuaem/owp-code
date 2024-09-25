%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uncertaintyAnalysis.m
% This script...
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 9/23/2024
% Last updated: /2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

%====Import data===========================================================
% Load results from duplicate sensor comparison
cd([rootpath,'open-water-platform-data\',site,'\cleaned\dupcheck'])
load([site,'-bestGuess.mat'])   % Output from dupCheck scripts

% Load results from Winkler validation
cd([rootpath,'open-water-platform-data\',site,'\cleaned\validation'])
load([site,'-winklerValid.mat'])
DOconc_mdl = mdl;

clear mdl

% Load results from DIC/TA discrete sample validation

% Load diel analysis results from different k parameterizations
cd([rootpath,'diel-method\sensitivity-analysis\',site,'\ro_hunt_wann_vary'])
R2b = load('diel_res.mat');
cd([rootpath,'diel-method\sensitivity-analysis\',site,'\wann_vary'])
W2b = load('diel_res.mat');
cd([rootpath,'diel-method\sensitivity-analysis\',site,'\emer_vary'])
Eb = load('diel_res.mat');

%====DO conc===============================================================
% Calculate standard error of the mean
SEM.DOconc = DOconc_mdl.RMSE / sqrt(height(DOconc_mdl.Variables));


%%
%====k parameterization====================================================
% Find average of daily standard deviations between k parameterizations
k_mean_std = mean(std([R2b.fas.k, W2b.fas.k],0,2),'omitmissing')

k_mean_std = mean(std([R2b.fas.k, Eb.fas.k],0,2),'omitmissing')

k_mean_std = mean(std([W2b.fas.k, Eb.fas.k],0,2),'omitmissing')