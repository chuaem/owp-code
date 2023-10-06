%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotAllDeployments.m
% This script plots the entire record of data (8 plots; one for each parameter)
% from one open-water platform for both sondes (BC and ERDC).
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% 9/27/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc

site = 'south'; % CHANGE THIS

cd(['G:\My Drive\Postdoc\SMIIL\raw-data\open-water-platform-data\',site])

ds = fileDatastore(['G:\My Drive\Postdoc\SMIIL\raw-data\open-water-platform-data\',site],"ReadFcn",@load,"FileExtensions",".mat");

dat = readall(ds);

sonde1_all = table();
sonde2_all = table();

for i = 1:length(dat)
    sonde1_all = [sonde1_all;dat{i}.sonde1];
end

for j = 1:length(dat)
    if strcmp(site,'south') == 1
        skipNum = [3,5];
        if ~ismember(j,skipNum)
            sonde2_all = [sonde2_all;dat{j}.sonde2];
        end
    else
        sonde2_all = [sonde2_all;dat{j}.sonde2];
    end
end

cd(['G:\My Drive\Postdoc\SMIIL\figures\open-water-platform-figures\',site])

% Global plotting settings
FontSize = 12;
NumTicks = 13;
LineWidth = 1;
XTick = linspace(sonde1_all.datetime_utc(1),sonde1_all.datetime_utc(end),NumTicks);
XTickFormat = "M/yy";
Legend = {'BC','ERDC'};
XLabel = 'Month/Year';
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde

figure(1),clf
plot(sonde1_all.datetime_utc,sonde1_all.depth,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.depth,'Color',blue)
hold off
ax1 = gcf().CurrentAxes;
ylabel('Depth (m)')
xlim('tight');ylim('tight')
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)

figure(2),clf
plot(sonde1_all.datetime_utc,sonde1_all.p,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.p,'Color',blue)
hold off
ylabel('Pressure (psi)')
xlim(ax1.XLim);ylim('tight')    % Use same x limits for all plots
ax2 = gcf().CurrentAxes;
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)

figure(3),clf
plot(sonde1_all.datetime_utc,sonde1_all.salinity,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.salinity,'Color',blue)
hold off
ylabel('Salinity (PSU)')
xlim(ax1.XLim);ylim('tight')    % Use same x limits for all plots
ax3 = gcf().CurrentAxes;
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)

figure(4),clf
plot(sonde1_all.datetime_utc,sonde1_all.temperature,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.temperature,'Color',blue)
hold off
ylabel('Temperature (^oC)')
xlim(ax1.XLim);ylim('tight')    % Use same x limits for all plots
ax4 = gcf().CurrentAxes;
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)

figure(5),clf
plot(sonde2_all.datetime_utc,sonde2_all.turbidity,'Color',blue)
ylabel('Turbidity (NTU)')
xlim(ax1.XLim);ylim('tight')    % Use same x limits for all plots
ax5 = gcf().CurrentAxes;
xtickformat(XTickFormat)
legend(Legend{2})
xlabel(XLabel)

figure(6),clf
plot(sonde1_all.datetime_utc,sonde1_all.pH,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.pH,'Color',blue)
hold off
ylabel('pH')
xlim(ax1.XLim);ylim('tight')    % Use same x limits for all plots
ax6 = gcf().CurrentAxes;
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)

figure(7),clf
plot(sonde1_all.datetime_utc,sonde1_all.DO_conc,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.DO_conc,'Color',blue)
hold off
ylabel('DO concentration (\mumol/L)')
xlim(ax1.XLim);ylim('tight')    % Use same x limits for all plots
ax7 = gcf().CurrentAxes;
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)

figure(8),clf
plot(sonde1_all.datetime_utc,sonde1_all.ORP,'Color',red)
hold on
plot(sonde2_all.datetime_utc,sonde2_all.ORP,'Color',blue)
hold off
ylabel('ORP (mV)')
xlim(ax1.XLim);ylim('tight')    % Use same x limits for all plots
ax8 = gcf().CurrentAxes;
xtickformat(XTickFormat)
legend(Legend)
xlabel(XLabel)

figure(9),clf
plot(sonde1_all.datetime_utc,sonde1_all.chla,'Color',red)
ylabel('Chl a (RFU)')
xlim(ax1.XLim);ylim('tight')    % Use same x limits for all plots
ax9 = gcf().CurrentAxes;
xtickformat(XTickFormat)
legend(Legend{1})
xlabel(XLabel)

% Set these properties for all figures
set([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9],'FontSize',FontSize,'LineWidth',LineWidth,'XTick',XTick)

cd(['G:\My Drive\Postdoc\SMIIL\figures\open-water-platform-figures\',site])