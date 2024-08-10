%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compareFinalQC.m
% This script compares the final QC'd results for all three sites.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 5/15/2024
% Last updated: 5/17/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

cd([rootpath,'open-water-platform-data\gull\cleaned\final-qc'])
load('gull-cleaned.mat');
gull = finalQC;

cd([rootpath,'open-water-platform-data\north\cleaned\final-qc'])
load('north-cleaned.mat');
north = finalQC;

cd([rootpath,'open-water-platform-data\south\cleaned\final-qc'])
load('south-cleaned.mat');
south = finalQC;

clearvars finalQC

% Find indices of deployment changes
ind_dep = find(diff(south.deployment) > 0);
dep = table([1;south.deployment(ind_dep+1)],[1;ind_dep+1]);
dep.Properties.VariableNames = {'depNum','ind'};

label = {'Deployment 1','Deployment 2','Deployment 4','Deployment 5',...
    'Deployment 6','Deployment 7','Deployment 8','Deployment 9',...
    'Deployment 10','Deployment 11','Deployment 12','Deployment 13',...
    'Deployment 14','Deployment 15','Deployment 16','Deployment 17','Deployment 18'};

fig1 = figure(1);clf
fig1.WindowState = 'maximized';
plot(gull.datetime_utc,gull.depth,'.','Color',rgb('goldenrod'),'DisplayName','Gull')
hold on
plot(north.datetime_utc,north.depth,'.','Color',rgb('lightseagreen'),'DisplayName','North')
plot(south.datetime_utc,south.depth,'.','Color',rgb('mediumpurple'),'DisplayName','South')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Depth (m)')
legend('show','location','best')
title('Final QC''d data')

fig2 = figure(2);clf
fig2.WindowState = 'maximized';
plot(gull.datetime_utc,gull.salinity,'.','Color',rgb('goldenrod'),'DisplayName','Gull')
hold on
plot(north.datetime_utc,north.salinity,'.','Color',rgb('lightseagreen'),'DisplayName','North')
plot(south.datetime_utc,south.salinity,'.','Color',rgb('mediumpurple'),'DisplayName','South')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Salinity (psu)')
legend('show','location','best')
title('Final QC''d data')

fig3 = figure(3);clf
fig3.WindowState = 'maximized';
plot(gull.datetime_utc,gull.temperature,'.','Color',rgb('goldenrod'),'DisplayName','Gull')
hold on
plot(north.datetime_utc,north.temperature,'.','Color',rgb('lightseagreen'),'DisplayName','North')
plot(south.datetime_utc,south.temperature,'.','Color',rgb('mediumpurple'),'DisplayName','South')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Temperature (^oC)')
legend('show','location','best')
title('Final QC''d data')

fig4 = figure(4);clf
fig4.WindowState = 'maximized';
plot(gull.datetime_utc,gull.DOconc,'.','Color',rgb('goldenrod'),'DisplayName','Gull')
hold on
plot(north.datetime_utc,north.DOconc,'.','Color',rgb('lightseagreen'),'DisplayName','North')
plot(south.datetime_utc,south.DOconc,'.','Color',rgb('mediumpurple'),'DisplayName','South')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('DO conc (\mumol/L)')
legend('show','location','best')
title('Final QC''d data')

fig5 = figure(5);clf
fig5.WindowState = 'maximized';
plot(gull.datetime_utc,gull.pH,'.','Color',rgb('goldenrod'),'DisplayName','Gull')
hold on
plot(north.datetime_utc,north.pH,'.','Color',rgb('lightseagreen'),'DisplayName','North')
plot(south.datetime_utc,south.pH,'.','Color',rgb('mediumpurple'),'DisplayName','South')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('pH')
legend('show','location','best')
title('Final QC''d data')

fig6 = figure(6);clf
fig6.WindowState = 'maximized';
plot(gull.datetime_utc,gull.turbidity,'.','Color',rgb('goldenrod'),'DisplayName','Gull')
hold on
plot(north.datetime_utc,north.turbidity,'.','Color',rgb('lightseagreen'),'DisplayName','North')
plot(south.datetime_utc,south.turbidity,'.','Color',rgb('mediumpurple'),'DisplayName','South')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Turbidity (NTU)')
legend('show','location','best')
title('Final QC''d data')

fig7 = figure(7);clf
fig7.WindowState = 'maximized';
plot(gull.datetime_utc,gull.chla,'.','Color',rgb('goldenrod'),'DisplayName','Gull')
hold on
plot(north.datetime_utc,north.chla,'.','Color',rgb('lightseagreen'),'DisplayName','North')
plot(south.datetime_utc,south.chla,'.','Color',rgb('mediumpurple'),'DisplayName','South')
xline(south.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
ylabel('Chl a (RFU)')
legend('show','location','best')
title('Final QC''d data')

%====Save the plots========================================================
option = questdlg('Save plots as .png and .fig?','Save plots','Yes','No','Yes');

switch option
    case 'Yes'
        cd([rootpath,'figures\open-water-platform\final-qc']);
    
        saveas(fig1,'finalQC_depth.png')
        saveas(fig1,'finalQC_depth.fig')
        saveas(fig2,'finalQC_salinity.png')
        saveas(fig2,'finalQC_salinity.fig')
        saveas(fig3,'finalQC_temperature.png')
        saveas(fig3,'finalQC_temperature.fig')
        saveas(fig4,'finalQC_DOconc.png')
        saveas(fig4,'finalQC_DOconc.fig')
        saveas(fig5,'finalQC_pH.png')
        saveas(fig5,'finalQC_pH.fig')
        saveas(fig6,'finalQC_turbidity.png')
        saveas(fig6,'finalQC_turbidity.fig')
        saveas(fig7,'finalQC_chla.png')
        saveas(fig7,'finalQC_chla.fig')

        disp('Plots saved!')

    case 'No'
        disp('Plots not saved.')
end