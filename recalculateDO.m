%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recalculateDO.m
% This script uses the "best-guess" salinity and temperature from
% dupCheck_d_S_T.m in the aaoptode_salpresscorr() function and recalculates 
% the DO concentration for each sonde.
%
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 4/26/2024
% Last updated: 5/15/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

%====Import sonde data=====================================================
% Best-guess data for depth, S, and T
cd([rootpath,'open-water-platform-data\',site,'\cleaned\dupcheck\'])
load([site,'-bestGuess.mat'])

% Other parameters
cd([rootpath,'open-water-platform-data\',site,'\cleaned\movmed\'])
load([site,'-bc-cleaned'])
load([site,'-erdc-cleaned'])
dat1 = sonde1_cleaned;
dat2 = sonde2_cleaned;
% Synchronize the BC and ERDC data to a common datetime vector
dat_syn = synchronize(dat1,dat2);

clearvars dat1 dat2 sonde1_cleaned sonde2_cleaned

% Find indices of deployment changes
ind_dep = find(diff(bestguess.deployment.deployment) > 0);
switch site
    case 'Gull'
        dep = table([1;dat_syn.deployment_dat1(ind_dep+1)],[1;ind_dep+1]);
    case 'North'
        dep = table([2;dat_syn.deployment_dat1(ind_dep+1)],[1;ind_dep+1]);
    case 'South'
        dep = table([1;dat_syn.deployment_dat1(ind_dep+1);15],[1;ind_dep+1;find(dat_syn.deployment_dat2 == 15,1)]);
        dep = sortrows(dep);
end
dep.Properties.VariableNames = {'depNum','ind'};

%====Recalculate DO========================================================
O2raw = dat_syn.DO_conc_dat1; % [uM]
temp = bestguess.temperature.temperature;   % [degC]
sal = bestguess.salinity.salinity;       % measured
press = bestguess.depth.depth;        % [dbar] (approx equal to [m])
S0 = dat_syn.salinity_dat1;   % internal salinity setting 
O2corr1 = aaoptode_salpresscorr(O2raw,temp,sal,press,S0); % [uM]

O2raw = dat_syn.DO_conc_dat2; % [uM]
temp = bestguess.temperature.temperature;   % [degC]
sal = bestguess.salinity.salinity;       % measured
press = bestguess.depth.depth;        % [dbar] (approx equal to [m])
S0 = dat_syn.salinity_dat2;   % internal salinity setting 
O2corr2 = aaoptode_salpresscorr(O2raw,temp,sal,press,S0); % [uM]

O2corr1 = table(dat_syn.datetime_utc,O2corr1);
O2corr1.Properties.VariableNames = {'datetime_utc','DOconc'};
O2corr1 = table2timetable(O2corr1);
O2corr1 = rmmissing(O2corr1);

O2corr2 = table(dat_syn.datetime_utc,O2corr2);
O2corr2.Properties.VariableNames = {'datetime_utc','DOconc'};
O2corr2 = table2timetable(O2corr2);
O2corr2 = rmmissing(O2corr2);

%====Global plotting settings==============================================
dt1 = datetime('29-Jun-2021','TimeZone','UTC');     % Make all plots have same start date
dt2 = dat_syn.datetime_utc(end);
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde
fontsize = 14;
linewidth = 1;
dotsize = 12;
circlesize = 6;

switch site
    case 'Gull'
        label = {'Deployment 1','Deployment 2','Deployment 5','Deployment 6',...
            'Deployment 7','Deployment 8','Deployment 9','Deployment 10',...
            'Deployment 11','Deployment 12','Deployment 13','Deployment 14',...
            'Deployment 15','Deployment 16','Deployment 17'};
    case 'North'
        label = {'Deployment 2','Deployment 6','Deployment 7','Deployment 8',...
            'Deployment 9','Deployment 10','Deployment 11','Deployment 12',...
            'Deployment 13','Deployment 14','Deployment 15','Deployment 16','Deployment 17'};
    case 'South'
        label = {'Deployment 1','Deployment 2','Deployment 4','Deployment 5',...
            'Deployment 6','Deployment 7','Deployment 8','Deployment 9',...
            'Deployment 10','Deployment 11','Deployment 12','Deployment 13',...
            'Deployment 14','Deployment 15','Deployment 16','Deployment 17'};
end

dt_utc = dat_syn.datetime_utc;

cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\DOconc'])
% Compare BC measured and recalculated DO conc
fig = figure(1);clf
fig.WindowState = 'maximized';
plot(dat_syn.datetime_utc,dat_syn.DO_conc_dat1,'.-','Color',red,'MarkerSize',dotsize,'DisplayName','BC - Measured')
hold on
plot(O2corr1.datetime_utc,O2corr1.DOconc,'.-','Color',rgb('darkred'),'MarkerSize',dotsize,'DisplayName','BC - Recalculated')
xline(dt_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('DO Concentration (\mumol/L)')
legend('show')
title([site,' BC - Measured and recalculated dissolved oxygen'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites

% Compare ERDC measured and recalculated DO conc
fig = figure(2);clf
fig.WindowState = 'maximized';
plot(dat_syn.datetime_utc,dat_syn.DO_conc_dat2,'.-','Color',blue,'MarkerSize',dotsize,'DisplayName','ERDC - Measured')
hold on
plot(O2corr2.datetime_utc,O2corr2.DOconc,'.-','Color',rgb('darkblue'),'MarkerSize',dotsize,'DisplayName','ERDC - Recalculated')
xline(dt_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('DO Concentration (\mumol/L)')
legend('show')
title([site,' ERDC - Measured and recalculated dissolved oxygen'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites

DOcorr = struct('bc',O2corr1,'erdc',O2corr2);

%====Save recalculated data===============================================
option = questdlg('Save recalculated data?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'open-water-platform-data\',site,'\cleaned\corrected'])
        save([site,'-DOcorr.mat'],'DOcorr')
        disp('File saved!')
    case 'No'
        disp('File not saved.')
end